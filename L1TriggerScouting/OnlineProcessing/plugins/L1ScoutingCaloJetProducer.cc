#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "CommonTools/Utils/interface/FormulaEvaluator.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloJet.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "EventFilter/L1ScoutingRawToDigi/interface/masks.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "L1TriggerScouting/Utilities/interface/conversion.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

class L1ScoutingCaloJetProducer : public edm::global::EDProducer<> {
public:
  explicit L1ScoutingCaloJetProducer(const edm::ParameterSet&);
  ~L1ScoutingCaloJetProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  class JetCorrector {
  public:
    explicit JetCorrector() = default;
    explicit JetCorrector(std::string const& filePath) {
      std::ifstream infile(filePath);
      if (!infile) {
        throw cms::Exception("InputError") << "failed to open JetCorrector input file: " << filePath;
      }

      std::string line{};
      while (std::getline(infile, line)) {
        std::istringstream iss(line);

        float ptMin{0.f};
        float ptMax{0.f};
        float etaMin{0.f};
        float etaMax{0.f};
        int puProxyMin{0};
        int puProxyMax{0};
        int formNParams{0};
        std::string formEvalStr{""};

        if (!(iss >> ptMin >> ptMax >> etaMin >> etaMax >> puProxyMin >> puProxyMax >> formEvalStr >> formNParams)) {
          throw cms::Exception("InvalidInput")
              << "failed to read line from input file (invalid format): \"" << line << "\"";
        }

        if (ptMin <= 0) {
          throw cms::Exception("InvalidInput")
              << "invalid value for parameter \"ptMin\" (must be greater than zero): " << ptMin;
        }

        if (ptMax <= 0) {
          throw cms::Exception("InvalidInput")
              << "invalid value for parameter \"ptMax\" (must be greater than zero): " << ptMax;
        }

        if (ptMin >= ptMax) {
          throw cms::Exception("InvalidInput") << "inconsistent values for parameters \"ptMin\" and \"ptMax\" (the "
                                                  "latter must be greater than the former): ptMin="
                                               << ptMin << " ptMax=" << ptMax;
        }

        if (etaMin >= etaMax) {
          throw cms::Exception("InvalidInput") << "inconsistent values for parameters \"etaMin\" and \"etaMax\" (the "
                                                  "latter must be greater than the former): etaMin="
                                               << etaMin << " etaMax=" << etaMax;
        }

        if (puProxyMin > 0 and puProxyMax > 0 and puProxyMin >= puProxyMax) {
          throw cms::Exception("InvalidInput")
              << "inconsistent values for parameters \"puProxyMin\" and \"puProxyMax\" (if both are greater than zero, "
                 "the latter must be greater than the former): puProxyMin="
              << puProxyMin << " puProxyMax=" << puProxyMax;
        }

        if (formNParams <= 0) {
          throw cms::Exception("InvalidInput")
              << "invalid value for parameter \"formNParams\" (must be greater than zero): " << formNParams;
        }

        reco::FormulaEvaluator formEval{formEvalStr};

        std::vector<double> formParams(formNParams);
        for (auto idx = 0; idx < formNParams; ++idx) {
          if (!(iss >> formParams[idx])) {
            throw cms::Exception("InvalidInput")
                << "failed to read line from input file (invalid format, formula parameter #" << idx << "): \"" << line
                << "\"";
          }
        }

        data_.emplace_back(
            ptMin, ptMax, etaMin, etaMax, puProxyMin, puProxyMax, std::move(formEval), std::move(formParams));
      }
    }

    double correction(float const pt, float const eta, int const puProxy) const {
      for (auto const& entry : data_) {
        if (eta >= entry.etaMin and eta < entry.etaMax and (entry.puProxyMin < 0 or puProxy >= entry.puProxyMin) and
            (entry.puProxyMax < 0 or puProxy < entry.puProxyMax)) {
          std::vector<double> vars{std::clamp(pt, entry.ptMin, entry.ptMax)};
          return (entry.formulaEvaluator.evaluate(vars, entry.formulaParameters) / vars[0]);
        }
      }
      return 0;
    }

  private:
    struct Entry {
      float ptMin;
      float ptMax;
      float etaMin;
      float etaMax;
      int puProxyMin;
      int puProxyMax;
      reco::FormulaEvaluator formulaEvaluator;
      std::vector<double> formulaParameters;
    };

    std::vector<Entry> data_;
  };

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // Number of BXs per orbit
  static constexpr unsigned int kNBXPlus1 = 3565;

  edm::EDGetTokenT<l1ScoutingRun3::CaloTowerOrbitCollection> const src_;
  double const akR_;
  double const ptMin_;
  int const towerMinHwEt_;
  int const towerMaxHwEt_;

  bool const applyJECs_;
  JetCorrector const jetCorrector_;
  int const jecPUProxyTowerMinHwEt_;
  int const jecPUProxyTowerMaxHwEt_;
  int const jecPUProxyTowerMinAbsHwEta_;
  int const jecPUProxyTowerMaxAbsHwEta_;

  bool const produceSortedCaloTowers_;

  int const mantissaPrecision_;
};

L1ScoutingCaloJetProducer::L1ScoutingCaloJetProducer(const edm::ParameterSet& iPSet)
    : src_(consumes(iPSet.getParameter<edm::InputTag>("src"))),
      akR_(iPSet.getParameter<double>("akR")),
      ptMin_(iPSet.getParameter<double>("ptMin")),
      towerMinHwEt_(iPSet.getParameter<int>("towerMinHwEt")),
      towerMaxHwEt_(iPSet.getParameter<int>("towerMaxHwEt")),
      applyJECs_(iPSet.getParameter<bool>("applyJECs")),
      jetCorrector_(applyJECs_ ? JetCorrector(iPSet.getParameter<edm::FileInPath>("jecFile").fullPath())
                               : JetCorrector{}),
      jecPUProxyTowerMinHwEt_(iPSet.getParameter<int>("jecPUProxyTowerMinHwEt")),
      jecPUProxyTowerMaxHwEt_(iPSet.getParameter<int>("jecPUProxyTowerMaxHwEt")),
      jecPUProxyTowerMinAbsHwEta_(iPSet.getParameter<int>("jecPUProxyTowerMinAbsHwEta")),
      jecPUProxyTowerMaxAbsHwEta_(iPSet.getParameter<int>("jecPUProxyTowerMaxAbsHwEta")),
      produceSortedCaloTowers_(iPSet.getParameter<bool>("produceSortedCaloTowers")),
      mantissaPrecision_(iPSet.getParameter<int>("mantissaPrecision")) {
  produces<l1ScoutingRun3::CaloJetOrbitCollection>("CaloJet").setBranchAlias("CaloJetOrbitCollection");
  if (produceSortedCaloTowers_) {
    produces<l1ScoutingRun3::CaloTowerOrbitCollection>("SortedCaloTowers");
  }
}

// ------------ method called for each ORBIT  ------------
void L1ScoutingCaloJetProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  auto const& caloTowerCollection = iEvent.get(src_);

  // Output containers for CaloJets
  auto caloJetCollection = std::make_unique<l1ScoutingRun3::CaloJetOrbitCollection>();
  std::vector<std::vector<l1ScoutingRun3::CaloJet>> caloJetBuffer(kNBXPlus1);
  unsigned int nCaloJet = 0;

  // Output containers for sorted CaloTowers (used only if "produceSortedCaloTowers == True")
  auto sortedCaloTowerCollection = std::make_unique<l1ScoutingRun3::CaloTowerOrbitCollection>();
  std::vector<std::vector<l1ScoutingRun3::CaloTower>> sortedCaloTowerBuffer(kNBXPlus1);
  unsigned int nSortedCaloTower = 0;

  // Define fastjet algorithm
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, akR_);

  // Loop over valid bunch crossings
  for (auto const bx : caloTowerCollection.getFilledBxs()) {
    LogTrace("L1ScoutingCaloJetProducer")
        << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel() << "] BX = " << bx;

    LogTrace("L1ScoutingCaloJetProducer") << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel()
                                          << "]   Inputs (l1ScoutingRun3::CaloTower and fastjet::PseudoJet)";

    auto const& cts = caloTowerCollection.bxIterator(bx);
    nSortedCaloTower += cts.size();

    // Indices of the CaloTowers that will not be assigned to any jets.
    // It starts as a vector containing all the values from 0 to "cts.size() - 1",
    // then the indices of CaloTowers used for jet clustering are removed one at a time.
    std::vector<int> unclusteredCaloTowerIndices(cts.size());
    std::iota(unclusteredCaloTowerIndices.begin(), unclusteredCaloTowerIndices.end(), 0);

    // Prepare PseudoJets to give in input to fastjet
    std::vector<fastjet::PseudoJet> pjCTs;
    pjCTs.reserve(cts.size());

    for (auto cidx{0u}; cidx < cts.size(); ++cidx) {
      auto const& ct{cts[cidx]};

      if (not((towerMinHwEt_ < 0 or ct.hwEt() >= towerMinHwEt_) and
              (towerMaxHwEt_ < 0 or ct.hwEt() <= towerMaxHwEt_))) {
        continue;
      }

      if (not l1ScoutingRun3::calol1::validHwEta(ct.hwEta())) {
        edm::LogWarning("L1ScoutingCaloJetProducer") << "CaloTower in BX=" << bx << " with invalid hwEta value ("
                                                     << ct.hwEta() << ") will not be used for jet clustering !";
        continue;
      }

      if (not l1ScoutingRun3::calol1::validHwPhi(ct.hwPhi())) {
        edm::LogWarning("L1ScoutingCaloJetProducer") << "CaloTower in BX=" << bx << " with invalid hwPhi value ("
                                                     << ct.hwPhi() << ") will not be used for jet clustering !";
        continue;
      }

      float const ctEt = l1ScoutingRun3::calol1::fEt(ct.hwEt());
      float const ctEta = l1ScoutingRun3::calol1::fEta(ct.hwEta());
      float const ctPhi = l1ScoutingRun3::calol1::fPhi(ct.hwPhi());

      pjCTs.emplace_back(fastjet::PtYPhiM(ctEt, ctEta, ctPhi, 0));
      pjCTs.back().set_user_index(cidx);

      LogTrace("L1ScoutingCaloJetProducer")
          << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel() << "]     [" << (pjCTs.size() - 1)
          << "] hwEt=" << ct.hwEt() << " hwEta=" << ct.hwEta() << " hwPhi=" << ct.hwPhi() << " (PseudoJet: pt=" << ctEt
          << " eta=" << ctEta << " phi=" << ctPhi << " px=" << pjCTs.back().px() << " py=" << pjCTs.back().py()
          << " pz=" << pjCTs.back().pz() << " E=" << pjCTs.back().E() << " user_index=" << pjCTs.back().user_index()
          << ")";
    }

    // If JECs are applied, compute a per-BX PU proxy used as input to the evaluation of the JECs
    // (the PU proxy corresponds to the number of CaloTowers passing predefined cuts on hwEt and |hwEta|)
    int puProxy{0};
    if (applyJECs_) {
      for (auto const& ct : cts) {
        if (not((jecPUProxyTowerMinHwEt_ < 0 or ct.hwEt() >= jecPUProxyTowerMinHwEt_) and
                (jecPUProxyTowerMaxHwEt_ < 0 or ct.hwEt() <= jecPUProxyTowerMaxHwEt_))) {
          continue;
        }

        auto const absHwEta{std::abs(ct.hwEta())};
        if (not((jecPUProxyTowerMinAbsHwEta_ < 0 or absHwEta >= jecPUProxyTowerMinAbsHwEta_) and
                (jecPUProxyTowerMaxAbsHwEta_ < 0 or absHwEta <= jecPUProxyTowerMaxAbsHwEta_))) {
          continue;
        }

        ++puProxy;
      }

      LogTrace("L1ScoutingCaloJetProducer") << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel()
                                            << "]   PU proxy for JECs (CaloTower multiplicity): " << puProxy;
    }

    LogTrace("L1ScoutingCaloJetProducer") << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel()
                                          << "]   Running jet clustering and applying JECs";

    // Run the jet clustering with the given jet definition
    fastjet::ClusterSequence clustSeq(pjCTs, jetDef);

    // Get the resulting jets ordered in pt
    std::vector<fastjet::PseudoJet> incJets = clustSeq.inclusive_jets();

    // Fill l1ScoutingRun3::CaloJet objects buffer
    std::vector<l1ScoutingRun3::CaloJet> unsortedCaloJets{};
    unsortedCaloJets.reserve(incJets.size());

    std::vector<std::vector<int>> unsortedCaloJetConstIndices{};
    unsortedCaloJetConstIndices.reserve(incJets.size());

    for (auto idx = 0u; idx < incJets.size(); ++idx) {
      auto const& incJet = incJets[idx];

      double const energyCorr{applyJECs_ ? jetCorrector_.correction(incJet.pt(), incJet.eta(), puProxy) : 1};

      LogTrace("L1ScoutingCaloJetProducer")
          << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel() << "]     [" << idx
          << "] pt=" << incJet.pt() << " eta=" << incJet.eta() << " phi=" << incJet.phi_std() << " mass=" << incJet.m()
          << " energyCorr=" << energyCorr << " (before JECs and pT cut)";

      if (energyCorr <= 0) {
        continue;
      }

      float const jet_pt = incJet.pt() * energyCorr;

      if (jet_pt < ptMin_) {
        continue;
      }

      float const jet_mass = incJet.m() * energyCorr;

      // Variables related to jet constituents
      int nConst{0};
      int nConstSatECAL{0};
      int nConstSatHCAL{0};
      int nConstSatECALAndHCAL{0};

      // "energyEm" ("energyTot") corresponds to the sum of the energies
      // measured in ECAL (ECAL + HCAL) of the CaloTowers assigned to the jet.
      // The ratio "energyEm/energyTot" will then be used as the EM fraction of the jet energy.
      // "energyTot" is used as denominator of this fraction
      // to guarantee that the fraction is a value between 0 and 1.
      // Further down in this block, a warning is emitted if "energyTot"
      // differs by more than 5% from the uncorrected energy of the jet.
      auto energyEm{0.f};
      auto energyTot{0.f};

      unsortedCaloJetConstIndices.emplace_back();
      auto& jetConstIndices{unsortedCaloJetConstIndices.back()};

      if (incJet.has_constituents()) {
        nConst = incJet.constituents().size();
        jetConstIndices.reserve(nConst);

        for (auto const& jet_const : incJet.constituents()) {
          auto const ct_idx{jet_const.user_index()};
          auto const& ct{cts[ct_idx]};

          jetConstIndices.emplace_back(ct_idx);
          unclusteredCaloTowerIndices.erase(
              std::remove(unclusteredCaloTowerIndices.begin(), unclusteredCaloTowerIndices.end(), ct_idx),
              unclusteredCaloTowerIndices.end());

          // CaloTower transverse energy (hardware value)
          auto const ctHwEt{ct.hwEt()};

          // Counters of CaloTowers with saturated energy in ECAL and/or HCAL
          if (ctHwEt == l1t::CaloTools::kSatEcal) {
            ++nConstSatECAL;
          } else if (ctHwEt == l1t::CaloTools::kSatHcal) {
            ++nConstSatHCAL;
          } else if (ctHwEt == l1t::CaloTools::kSatTower) {
            ++nConstSatECAL;
            ++nConstSatHCAL;
            ++nConstSatECALAndHCAL;
          }

          // Energy-ratio bits
          uint8_t const ctHwEtRatio = ct.erBits() & l1ScoutingRun3::calol1::masksCaloTowers::erBits;

          // Special bits for "zero flag" and "e over h"
          bool const ctZeroFlag = ct.miscBits() & 0b01;
          bool const ctEohrFlag = ct.miscBits() & 0b10;

          // CaloTower energy (physical value)
          float const ctEnergy{l1ScoutingRun3::calol1::fEt(ctHwEt) *
                               std::cosh(l1ScoutingRun3::calol1::fEta(ct.hwEta()))};

          // ctEFracEm: EM/ECAL fraction of the CaloTower's energy
          //  - The criteria below to determine the value of "ctEFracEm" from the "zero flag" and "eoh flag"
          //    of the CaloTower's "miscBits" are based on the implementation of the L1T CaloLayer1 emulation.
          //    https://github.com/cms-sw/cmssw/blob/e3685ee38b3c50d33912a6e6817dc468e4ca1812/L1Trigger/L1TCaloLayer1/src/UCTTower.cc#L65-L84
          auto ctEFracEm{0.f};
          if (ctZeroFlag) {
            ctEFracEm = ctEohrFlag ? 1.f : 0.f;
          } else {
            float const frac{1.f / (1.f + (1 << ctHwEtRatio))};
            float const antifrac{1.f - frac};
            ctEFracEm = ctEohrFlag ? antifrac : frac;
          }

          energyEm += ctEnergy * ctEFracEm;
          energyTot += ctEnergy;
        }
      }

      float const energyFracEm{energyTot > 0 ? energyEm / energyTot : 0.f};

      // Order indices of the jet's constituents by hwEt, and
      // append them to the vector of indices of the output CaloTowers
      if (produceSortedCaloTowers_) {
        std::stable_sort(jetConstIndices.begin(), jetConstIndices.end(), [&cts](auto const idx1, auto const idx2) {
          return cts[idx1].hwEt() > cts[idx2].hwEt();
        });
      }

      LogTrace("L1ScoutingCaloJetProducer")
          << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel() << "]     [" << idx << "] pt=" << jet_pt
          << " eta=" << incJet.eta() << " phi=" << incJet.phi_std() << " mass=" << jet_mass
          << " energyCorr=" << energyCorr << " energyFracEm=" << energyFracEm
          << " (Eem+Ehad)/E=" << energyTot / incJet.E() << " nConst=" << nConst << " nConstSatECAL=" << nConstSatECAL
          << " nConstSatHCAL=" << nConstSatHCAL << " nConstSatECALAndHCAL=" << nConstSatECALAndHCAL
          << " (after JECs and pT cut)";

      // Emit a warning if the denominator of the EM energy fraction
      // differs by more than 5% from the uncorrected energy of the jet
      if (std::abs(energyTot - incJet.E()) > 1.05f * incJet.E()) {
        edm::LogWarning("L1ScoutingCaloJetProducer")
            << "sum of estimated ECAL+HCAL CaloTowers' energies (" << energyTot
            << ") differs from total uncorrected jet energy (" << incJet.E() << ") by more than 5%!";
      }

      unsortedCaloJets.emplace_back(
          MiniFloatConverter::reduceMantissaToNbitsRounding(jet_pt, mantissaPrecision_),
          MiniFloatConverter::reduceMantissaToNbitsRounding(incJet.eta(), mantissaPrecision_),
          MiniFloatConverter::reduceMantissaToNbitsRounding(incJet.phi_std(), mantissaPrecision_),
          MiniFloatConverter::reduceMantissaToNbitsRounding(jet_mass, mantissaPrecision_),
          MiniFloatConverter::reduceMantissaToNbitsRounding(energyCorr, mantissaPrecision_),
          MiniFloatConverter::reduceMantissaToNbitsRounding(energyFracEm, mantissaPrecision_),
          nConst,
          nConstSatECAL,
          nConstSatHCAL,
          nConstSatECALAndHCAL);
      ++nCaloJet;
    }

    std::vector<int> caloJetIndices(unsortedCaloJets.size());
    std::iota(caloJetIndices.begin(), caloJetIndices.end(), 0);

    std::stable_sort(
        caloJetIndices.begin(), caloJetIndices.end(), [&unsortedCaloJets](auto const idx1, auto const idx2) {
          return unsortedCaloJets[idx1].pt() > unsortedCaloJets[idx2].pt();
        });

    LogTrace("L1ScoutingCaloJetProducer") << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel()
                                          << "]   Final outputs (l1ScoutingRun3::CaloJet)";

    auto& bxCaloJetBuffer = caloJetBuffer[bx];
    bxCaloJetBuffer.reserve(caloJetIndices.size());

    for (auto const idx : caloJetIndices) {
      bxCaloJetBuffer.emplace_back(unsortedCaloJets[idx]);
    }

    if (produceSortedCaloTowers_) {
      auto& bxSortedCaloTowerBuffer = sortedCaloTowerBuffer[bx];
      bxSortedCaloTowerBuffer.reserve(cts.size());

      for (auto const idx : caloJetIndices) {
        for (auto const idx2 : unsortedCaloJetConstIndices[idx]) {
          bxSortedCaloTowerBuffer.emplace_back(cts[idx2]);
        }
      }

      for (auto const idx2 : unclusteredCaloTowerIndices) {
        bxSortedCaloTowerBuffer.emplace_back(cts[idx2]);
      }
    }

#ifdef EDM_ML_DEBUG
    for (auto idx0{0u}; idx0 < caloJetIndices.size(); ++idx0) {
      auto const idx = caloJetIndices[idx0];
      auto const& obj = unsortedCaloJets[idx];
      LogTrace("L1ScoutingCaloJetProducer")
          << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel() << "]     [" << idx0
          << "] index=" << idx << " pt=" << obj.pt() << " eta=" << obj.eta() << " phi=" << obj.phi()
          << " mass=" << obj.mass() << " energyCorr=" << obj.energyCorr() << " energyFracEm=" << obj.energyFracEm()
          << " nConst=" << obj.nConst() << " nConstSatEnergyECAL=" << obj.nConstSaturatedEnergyECAL()
          << " nConstSatEnergyHCAL=" << obj.nConstSaturatedEnergyHCAL()
          << " nConstSatEnergyECALAndHCAL=" << obj.nConstSaturatedEnergyECALAndHCAL();

      if (produceSortedCaloTowers_) {
        for (auto const idx2 : unsortedCaloJetConstIndices[idx]) {
          auto const& ct{cts[idx2]};
          LogTrace("L1ScoutingCaloJetProducer")
              << "[L1ScoutingCaloJetProducer:" << moduleDescription().moduleLabel() << "]          CaloTower[" << idx2
              << "] hwEt=" << ct.hwEt() << " hwEta=" << ct.hwEta() << " hwPhi=" << ct.hwPhi()
              << " erBits=" << ct.erBits() << " miscBits=" << ct.miscBits();
        }
      }
    }
#endif
  }

  // fill orbit collection with reconstructed jets
  caloJetCollection->fillAndClear(caloJetBuffer, nCaloJet);
  iEvent.put(std::move(caloJetCollection), "CaloJet");

  if (produceSortedCaloTowers_) {
    sortedCaloTowerCollection->fillAndClear(sortedCaloTowerBuffer, nSortedCaloTower);
    iEvent.put(std::move(sortedCaloTowerCollection), "SortedCaloTowers");
  }
}

void L1ScoutingCaloJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src")->setComment(
      "Input collection of CaloTowers (type: l1ScoutingRun3::CaloTowerOrbitCollection)");
  desc.add<double>("akR")->setComment("Value of R parameter for anti-kt clustering");
  desc.add<double>("ptMin")->setComment(
      "Minimum pT of output jets (applies to the corrected jet-pT if applyJECs==True)");
  desc.add<int>("towerMinHwEt", 1)
      ->setComment("Min hwEt (inclusive) of CaloTowers used for jet clustering (ignored if negative)");
  desc.add<int>("towerMaxHwEt", -1)
      ->setComment("Max hwEt (inclusive) of CaloTowers used for jet clustering (ignored if negative)");

  desc.add<bool>("applyJECs", false)
      ->setComment("Apply jet-energy-scale corrections (and output corrected jets, ordered by their corrected pT)");
  desc.add<edm::FileInPath>("jecFile")->setComment(
      "Path to text file containing jet-energy-scale corrections (used only if applyJECs==True)");
  desc.add<int>("jecPUProxyTowerMinHwEt", 1)
      ->setComment(
          "Min CaloTower hwEt (inclusive) used when computing the CaloTower multiplicity taken as PU proxy to evaluate "
          "JECs (used only if applyJECs==True, and ignored if negative)");
  desc.add<int>("jecPUProxyTowerMaxHwEt", -1)
      ->setComment(
          "Max CaloTower hwEt (inclusive) used when computing the CaloTower multiplicity taken as PU proxy to evaluate "
          "JECs (used only if applyJECs==True, and ignored if negative)");
  desc.add<int>("jecPUProxyTowerMinAbsHwEta", 0)
      ->setComment(
          "Min CaloTower |hwEta| (inclusive) used when computing the CaloTower multiplicity taken as PU proxy to "
          "evaluate JECs (used only if applyJECs==True, and ignored if negative)");
  desc.add<int>("jecPUProxyTowerMaxAbsHwEta", 4)
      ->setComment(
          "Max CaloTower |hwEta| (inclusive) used when computing the CaloTower multiplicity taken as PU proxy to "
          "evaluate JECs (used only if applyJECs==True, and ignored if negative)");

  desc.add<bool>("produceSortedCaloTowers", false)
      ->setComment("Output a copy of the l1ScoutingRun3::CaloTowerOrbitCollection in \"src\" with a custom sorting");

  desc.add<int>("mantissaPrecision", 10)->setComment("default float16, change to 23 for float32");

  descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ScoutingCaloJetProducer);

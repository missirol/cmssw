#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

#define LOG edm::LogPrint("L1TCaloParamsViewer2") << "[L1TCaloParamsViewer2] "

class L1TCaloParamsViewer2 : public edm::one::EDAnalyzer<> {
public:
  explicit L1TCaloParamsViewer2(edm::ParameterSet const&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(edm::Event const&, edm::EventSetup const& iSetup) override;

  void printL1TCaloParamsHelper(l1t::CaloParamsHelper const& cph) const;

  edm::ESGetToken<l1t::CaloParams, L1TCaloParamsRcd> const caloParamsToken_;
  edm::ESWatcher<L1TCaloParamsRcd> caloParamsWatcher_;
};

L1TCaloParamsViewer2::L1TCaloParamsViewer2(edm::ParameterSet const& iConfig)
    : caloParamsToken_{esConsumes()}, caloParamsWatcher_{} {}

void L1TCaloParamsViewer2::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  if (caloParamsWatcher_.check(iSetup)) {
    printL1TCaloParamsHelper(iSetup.getData(caloParamsToken_));
  }
}

void L1TCaloParamsViewer2::printL1TCaloParamsHelper(l1t::CaloParamsHelper const& cph) const {
  LOG << "==========================================";
  LOG << "isValidForStage2 = " << cph.isValidForStage2();
  LOG << "regionLsb = " << cph.regionLsb();
  LOG << "etSumLsb = " << cph.etSumLsb();

  for (auto idx = 0u; idx < cph.etSumEtaMinSize(); ++idx) {
    LOG << "etSumEtaMin[" << idx << "] = " << cph.etSumEtaMin(idx);
  }

  for (auto idx = 0u; idx < cph.etSumEtaMaxSize(); ++idx) {
    LOG << "etSumEtaMax[" << idx << "] = " << cph.etSumEtaMax(idx);
  }

  for (auto idx = 0u; idx < cph.etSumEtThresholdSize(); ++idx) {
    LOG << "etSumEtThreshold[" << idx << "] = " << cph.etSumEtThreshold(idx);
  }

  LOG << "---------------------";
  LOG << "towerLsbH = " << cph.towerLsbH();
  LOG << "towerLsbE = " << cph.towerLsbE();
  LOG << "towerLsbSum = " << cph.towerLsbSum();
  LOG << "towerNBitsH = " << cph.towerNBitsH();
  LOG << "towerNBitsE = " << cph.towerNBitsE();
  LOG << "towerNBitsSum = " << cph.towerNBitsSum();
  LOG << "towerNBitsRatio = " << cph.towerNBitsRatio();
  LOG << "towerMaskE = " << cph.towerMaskE();
  LOG << "towerMaskH = " << cph.towerMaskH();
  LOG << "towerMaskSum = " << cph.towerMaskSum();
  LOG << "towerMaskRatio = " << cph.towerMaskRatio();
  LOG << "doTowerEncoding = " << cph.doTowerEncoding();

  LOG << "---------------------";
  LOG << "egLsb = " << cph.egLsb();
  LOG << "egSeedThreshold = " << cph.egSeedThreshold();
  LOG << "egNeighbourThreshold = " << cph.egNeighbourThreshold();
  LOG << "egHcalThreshold = " << cph.egHcalThreshold();
  LOG << "egMaxHcalEt = " << cph.egMaxHcalEt();
  LOG << "egMaxPtHOverE = " << cph.egMaxPtHOverE();
  LOG << "egMinPtJetIsolation = " << cph.egMinPtJetIsolation();
  LOG << "egMaxPtJetIsolation = " << cph.egMaxPtJetIsolation();
  LOG << "egMinPtHOverEIsolation = " << cph.egMinPtHOverEIsolation();
  LOG << "egMaxPtHOverEIsolation = " << cph.egMaxPtHOverEIsolation();
  LOG << "egIsoAreaNrTowersEta = " << cph.egIsoAreaNrTowersEta();
  LOG << "egIsoAreaNrTowersPhi = " << cph.egIsoAreaNrTowersPhi();
  LOG << "egIsoVetoNrTowersPhi = " << cph.egIsoVetoNrTowersPhi();

  LOG << "---------------------";
  LOG << "tauLsb = " << cph.tauLsb();
  LOG << "tauSeedThreshold = " << cph.tauSeedThreshold();
  LOG << "tauNeighbourThreshold = " << cph.tauNeighbourThreshold();
  LOG << "tauMaxPtTauVeto = " << cph.tauMaxPtTauVeto();
  LOG << "tauMinPtJetIsolationB = " << cph.tauMinPtJetIsolationB();
  LOG << "tauMaxJetIsolationB = " << cph.tauMaxJetIsolationB();
  LOG << "tauMaxJetIsolationA = " << cph.tauMaxJetIsolationA();
  LOG << "isoTauEtaMin = " << cph.isoTauEtaMin();
  LOG << "isoTauEtaMax = " << cph.isoTauEtaMax();
  LOG << "tauIsoAreaNrTowersEta = " << cph.tauIsoAreaNrTowersEta();
  LOG << "tauIsoAreaNrTowersPhi = " << cph.tauIsoAreaNrTowersPhi();
  LOG << "tauIsoVetoNrTowersPhi = " << cph.tauIsoVetoNrTowersPhi();

  LOG << "---------------------";
  LOG << "jetLsb = " << cph.jetLsb();
  LOG << "jetSeedThreshold = " << cph.jetSeedThreshold();
  LOG << "jetNeighbourThreshold = " << cph.jetNeighbourThreshold();

  auto const& nodes = cph.getNodes();
  unsigned int node_idx{0};
  for (auto const& node : nodes) {
    LOG << "---------------------";
    LOG << "Node " << node_idx;
    LOG << "  type = " << node.type_;
    LOG << "  version = " << node.version_;
    LOG << "  LUT.nrBitsAddress = " << node.LUT_.nrBitsAddress();
    LOG << "  LUT.nrBitsData = " << node.LUT_.nrBitsData();
    LOG << "  LUT.maxSize = " << node.LUT_.maxSize();
    LOG << "  LUT.empty = " << node.LUT_.empty();
    for (auto idx = 0u; idx < node.LUT_.maxSize(); ++idx) {
      LOG << "  LUT.data(" << idx << ") = " << node.LUT_.data(idx);
    }
    for (auto idx = 0u; idx < node.dparams_.size(); ++idx) {
      LOG << "  dparams[" << idx << "] = " << node.dparams_[idx];
    }
    for (auto idx = 0u; idx < node.uparams_.size(); ++idx) {
      LOG << "  uparams[" << idx << "] = " << node.uparams_[idx];
    }
    for (auto idx = 0u; idx < node.iparams_.size(); ++idx) {
      LOG << "  iparams[" << idx << "] = " << node.iparams_[idx];
    }
    for (auto idx = 0u; idx < node.sparams_.size(); ++idx) {
      LOG << "  sparams[" << idx << "] = " << node.sparams_[idx];
    }
    ++node_idx;
  }

  LOG << "==========================================";
}

void L1TCaloParamsViewer2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  descriptions.add("l1tCaloParamsViewer2", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TCaloParamsViewer2);

#ifndef L1TriggerScouting_OnlineProcessing_L1ScoutingCaloJetClusterizer_h
#define L1TriggerScouting_OnlineProcessing_L1ScoutingCaloJetClusterizer_h

#include <span>
#include <string>
#include <vector>

#include "CommonTools/Utils/interface/FormulaEvaluator.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloJet.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCaloTower.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "fastjet/JetDefinition.hh"

class L1ScoutingCaloJetClusterizer {
public:
  explicit L1ScoutingCaloJetClusterizer(const edm::ParameterSet&);

  static void fillDescription(edm::ParameterSetDescription& desc);

  void run(std::span<const l1ScoutingRun3::CaloTower> const cts,
           std::vector<l1ScoutingRun3::CaloJet>& bxCaloJetBuffer,
           std::vector<l1ScoutingRun3::CaloTower>& bxSortedCaloTowerBuffer,
           bool const fillSortedCaloTowers) const;

private:
  class JetCorrector {
  public:
    explicit JetCorrector() = default;
    explicit JetCorrector(std::string const& filePath);

    double correction(float const pt, float const eta, int const puProxy) const;

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

  fastjet::JetDefinition const jetDef_;
  double const ptMin_;
  int const towerMinHwEt_;
  int const towerMaxHwEt_;

  bool const applyJECs_;
  JetCorrector const jetCorrector_;

  int const jecPUProxyTowerMinHwEt_;
  int const jecPUProxyTowerMaxHwEt_;
  int const jecPUProxyTowerMinAbsHwEta_;
  int const jecPUProxyTowerMaxAbsHwEta_;

  int const mantissaPrecision_;
};

#endif

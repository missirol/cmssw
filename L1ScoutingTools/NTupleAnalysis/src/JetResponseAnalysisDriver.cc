#include <cmath>

#include "L1ScoutingTools/NTupleAnalysis/interface/JetResponseAnalysisDriver.h"
#include "L1ScoutingTools/NTupleAnalysis/interface/Utils.h"

void JetResponseAnalysisDriver::init() {
  // hard-coded for now..
  jecA_.init("/eos/cms/store/cmst3/group/daql1scout/run3_calotowers/jet_pt_corrections/mc_qcd_2025/graph_SC.root");

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT000to040"] = [](float a, int b) { return a <= 0.2 and b < 40; };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT000to040"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT000to040"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT000to040"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT000to040"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT000to040"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT000to040"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT000to040"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT000to040"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and b < 40;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT000to040"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and b < 40;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT040to080"] = [](float a, int b) {
    return a <= 0.2 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT040to080"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT040to080"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT040to080"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT040to080"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT040to080"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT040to080"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT040to080"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT040to080"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 40 <= b and b < 80;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT040to080"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 40 <= b and b < 80;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT080to120"] = [](float a, int b) {
    return a <= 0.2 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT080to120"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT080to120"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT080to120"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT080to120"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT080to120"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT080to120"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT080to120"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT080to120"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 80 <= b and b < 120;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT080to120"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 80 <= b and b < 120;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT120to160"] = [](float a, int b) {
    return a <= 0.2 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT120to160"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT120to160"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT120to160"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT120to160"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT120to160"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT120to160"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT120to160"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT120to160"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 120 <= b and b < 160;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT120to160"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 120 <= b and b < 160;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT160to200"] = [](float a, int b) {
    return a <= 0.2 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT160to200"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT160to200"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT160to200"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT160to200"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT160to200"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT160to200"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT160to200"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT160to200"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 160 <= b and b < 200;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT160to200"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 160 <= b and b < 200;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT200to240"] = [](float a, int b) {
    return a <= 0.2 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT200to240"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT200to240"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT200to240"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT200to240"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT200to240"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT200to240"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT200to240"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT200to240"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 200 <= b and b < 240;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT200to240"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 200 <= b and b < 240;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT240to300"] = [](float a, int b) {
    return a <= 0.2 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT240to300"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT240to300"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT240to300"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT240to300"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT240to300"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT240to300"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT240to300"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT240to300"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 240 <= b and b < 300;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT240to300"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 240 <= b and b < 300;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT300to400"] = [](float a, int b) {
    return a <= 0.2 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT300to400"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT300to400"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT300to400"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT300to400"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT300to400"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT300to400"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT300to400"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT300to400"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 300 <= b and b < 400;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT300to400"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 300 <= b and b < 400;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT400to600"] = [](float a, int b) {
    return a <= 0.2 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT400to600"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT400to600"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT400to600"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT400to600"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT400to600"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT400to600"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT400to600"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT400to600"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 400 <= b and b < 600;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT400to600"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 400 <= b and b < 600;
  };

  jetCategoryForJECFuncMap_["_absEta0p0to0p2nCT600toInf"] = [](float a, int b) { return a <= 0.2 and 600 <= b; };
  jetCategoryForJECFuncMap_["_absEta0p2to0p4nCT600toInf"] = [](float a, int b) {
    return 0.2 < a and a <= 0.4 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta0p4to0p6nCT600toInf"] = [](float a, int b) {
    return 0.4 < a and a <= 0.6 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta0p6to0p8nCT600toInf"] = [](float a, int b) {
    return 0.6 < a and a <= 0.8 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta0p8to1p0nCT600toInf"] = [](float a, int b) {
    return 0.8 < a and a <= 1.0 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta1p0to1p3nCT600toInf"] = [](float a, int b) {
    return 1.0 < a and a <= 1.3 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta1p3to1p6nCT600toInf"] = [](float a, int b) {
    return 1.3 < a and a <= 1.6 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta1p6to1p9nCT600toInf"] = [](float a, int b) {
    return 1.6 < a and a <= 1.9 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta1p9to2p5nCT600toInf"] = [](float a, int b) {
    return 1.9 < a and a <= 2.5 and 600 <= b;
  };
  jetCategoryForJECFuncMap_["_absEta2p5to3p0nCT600toInf"] = [](float a, int b) {
    return 2.5 < a and a <= 3.0 and 600 <= b;
  };

  // histogram: events counter
  addTH1D("eventsProcessed", {0, 1});
  addTH1D("weight", 100, -5, 5);
  addTH1D("nPU", 40, 0, 120);
  addTH1D("nCT", 100, 0, 1000);
  addTH2D("nPU__vs__nCT", 40, 0, 120, 100, 0, 1000);

  labelMap_jetAK4_ = {
    {"L1EmulJet", {{"GEN", "GenJet"}}},
    {"L1EmulAK4CTJet0", {{"GEN", "GenJet"}}},
    {"L1EmulAK4CTJet0CorrA", {{"GEN", "GenJet"}}},
  };

  for (auto const& selLabel : {"NoSelection"}) {
    // histograms: AK4 Jets
    for (auto const& jetLabel : labelMap_jetAK4_) {
      bookHistograms_Jets(selLabel, jetLabel.first, utils::mapKeys(jetLabel.second));
    }
  }
}

std::vector<std::string> JetResponseAnalysisDriver::jetCategoryLabelsForJECHistos(
    const std::string& jetColl, const std::string& matchJetLabel) const {
  std::vector<std::string> ret;
  ret.reserve(jetCategoryForJECFuncMap_.size());
  for (auto const& [key, foo] : jetCategoryForJECFuncMap_) {
    ret.emplace_back(key);
  }

  return ret;
}

void JetResponseAnalysisDriver::analyze() {
  H1("eventsProcessed")->Fill(0.5);

  float const wgt{1.f};
  H1("weight")->Fill(wgt);

  auto const nPU = this->value<float>("Pileup_nTrueInt");
  H1("nPU")->Fill(nPU, wgt);

  auto const nCT = this->value<int>("nL1EmulCaloTower");
  H1("nCT")->Fill(nCT, wgt);

  H2("nPU__vs__nCT")->Fill(nPU, nCT, wgt);

  // AK4 Jets
  const float minAK4JetPt(10.);
  const float minAK4JetPtRef(5.);
  const float maxAK4JetDeltaRmatchRef(0.2);

  for (auto const& jetLabel : labelMap_jetAK4_) {
    fillHistoDataJets fhDataAK4Jets;
    fhDataAK4Jets.jetCollection = jetLabel.first;
    fhDataAK4Jets.jetPtMin = (jetLabel.first == "GenJet") ? minAK4JetPtRef : minAK4JetPt;
    fhDataAK4Jets.jetPtMax = (jetLabel.first == "L1EmulJet") ? 1023.4 : -1;
    fhDataAK4Jets.jetAbsEtaMax = 5.0;
    for (auto const& jetLabelRefs : jetLabel.second) {
      auto const jetPtMin2 = (jetLabelRefs.second == "GenJet") ? minAK4JetPtRef : minAK4JetPt;
      auto const jetPtMax2 = (jetLabelRefs.second == "L1EmulJet") ? 1023.4 : -1;
      fhDataAK4Jets.matches.emplace_back(fillHistoDataJets::Match(
          jetLabelRefs.first, jetLabelRefs.second, jetPtMin2, jetPtMax2, maxAK4JetDeltaRmatchRef));
    }

    fillHistograms_Jets("NoSelection", fhDataAK4Jets, wgt);
  }
}

void JetResponseAnalysisDriver::bookHistograms_Jets(const std::string& dir,
                                                          const std::string& jetType,
                                                          const std::vector<std::string>& matchLabels) {
  auto dirPrefix(dir);
  while (dirPrefix.back() == '/') {
    dirPrefix.pop_back();
  }
  if (not dirPrefix.empty()) {
    dirPrefix += "/";
  }

  std::vector<float> binEdges_pt(104);
  for (uint idx = 0; idx < binEdges_pt.size(); ++idx) {
    binEdges_pt.at(idx) = 10. * idx;
  }

  std::vector<float> binEdges_eta(101);
  for (uint idx = 0; idx < binEdges_eta.size(); ++idx) {
    binEdges_eta.at(idx) = -5.0 + 0.1 * idx;
  }

  std::vector<float> binEdges_response(101);
  for (uint idx = 0; idx < binEdges_response.size(); ++idx) {
    binEdges_response.at(idx) = 0.05 * idx;
  }

  std::vector<float> binEdges_nPU(41);
  for (uint idx = 0; idx < binEdges_nPU.size(); ++idx) {
    binEdges_nPU.at(idx) = 3 * idx;
  }

  std::vector<float> binEdges_nCT(101);
  for (uint idx = 0; idx < binEdges_nCT.size(); ++idx) {
    binEdges_nCT.at(idx) = 10 * idx;
  }

  for (auto const& matchLabel : matchLabels) {
    auto const jetCategoryLabelsForJECHistos_v = jetCategoryLabelsForJECHistos(jetType, matchLabel);
    for (auto const& catLabel : jetCategoryLabelsForJECHistos_v) {
      addTH2D(dirPrefix + jetType + catLabel + "_MatchedTo" + matchLabel + "_pt_over" + matchLabel + "__vs__pt", binEdges_response, binEdges_pt);
      addTH2D(dirPrefix + jetType + catLabel + "_MatchedTo" + matchLabel + "_pt_over" + matchLabel + "__vs__" + matchLabel + "_pt", binEdges_response, binEdges_pt);
      addTH1D(dirPrefix + jetType + catLabel + "_MatchedTo" + matchLabel + "_eta", binEdges_eta);
      addTH1D(dirPrefix + jetType + catLabel + "_MatchedTo" + matchLabel + "_nCT", binEdges_nCT);
      addTH1D(dirPrefix + jetType + catLabel + "_MatchedTo" + matchLabel + "_nPU", binEdges_nPU);
    }
  }
}

void JetResponseAnalysisDriver::fillHistograms_Jets(const std::string& dir, const fillHistoDataJets& fhData, float const weight) {
  auto dirPrefix(dir);
  while (dirPrefix.back() == '/') {
    dirPrefix.pop_back();
  }
  if (not dirPrefix.empty()) {
    dirPrefix += "/";
  }

  auto const jetCollRequiresJecA{utils::stringEndsWith(fhData.jetCollection, "CorrA")};
  auto const jetCollBranchName{jetCollRequiresJecA ? fhData.jetCollection.substr(0, fhData.jetCollection.size() - 5)
                                                   : fhData.jetCollection};

  auto const nPU = this->value<float>("Pileup_nTrueInt");
  auto const nCT = this->value<int>("nL1EmulCaloTower");

  if (not hasTTreeReaderValue("n" + jetCollBranchName)) {
    return;
  }

  auto const v_pt_size = this->value<int>("n" + jetCollBranchName);

  std::vector<float> v_pt{};
  std::vector<float> v_eta{};
  std::vector<float> v_phi{};

  v_pt.reserve(v_pt_size);
  v_eta.reserve(v_pt_size);
  v_phi.reserve(v_pt_size);

  auto const& a_pt = this->array<float>(jetCollBranchName + "_pt");
  auto const& a_eta = this->array<float>(jetCollBranchName + "_eta");
  auto const& a_phi = this->array<float>(jetCollBranchName + "_phi");

  for (auto idx = 0; idx < v_pt_size; ++idx) {
    float corr = 1;
    if (jetCollRequiresJecA) {
      corr = jecA_.correction(a_pt[idx], a_eta[idx]);
    }

    v_pt.emplace_back(a_pt[idx] * corr);
    v_eta.emplace_back(a_eta[idx]);
    v_phi.emplace_back(a_phi[idx]);
  }

  std::vector<size_t> fhDataIndices{};
  fhDataIndices.reserve(v_pt_size);
  for (auto idx = 0; idx < v_pt_size; ++idx) {
    auto const passesMinPtCut = fhData.jetPtMin < 0 or v_pt[idx] > fhData.jetPtMin;
    auto const passesMaxPtCut = fhData.jetPtMax < 0 or v_pt[idx] < fhData.jetPtMax;
    auto const passesAbsEtaCut = fhData.jetAbsEtaMax < 0 or std::abs(v_eta[idx]) < fhData.jetAbsEtaMax;
    if (passesMinPtCut and passesMaxPtCut and passesAbsEtaCut) {
      fhDataIndices.emplace_back(idx);
    }
  }

  for (auto const& fhDataMatch : fhData.matches) {
    auto const matchLabel(fhDataMatch.label);
    auto const matchJetColl(fhDataMatch.jetCollection);
    auto const matchJetPtMin(fhDataMatch.jetPtMin);
    auto const matchJetPtMax(fhDataMatch.jetPtMax);
    auto const matchJetDeltaR2Min{fhDataMatch.jetDeltaRMin * fhDataMatch.jetDeltaRMin};

    auto const matchJetCollRequiresJecA{utils::stringEndsWith(matchJetColl, "CorrA")};
    auto const matchJetCollBranchName{matchJetCollRequiresJecA ? matchJetColl.substr(0, matchJetColl.size() - 5)
                                                               : matchJetColl};

    if (not hasTTreeReaderValue("n" + matchJetCollBranchName)) {
      continue;
    }

    auto const v_match_pt_size = this->value<int>("n" + matchJetCollBranchName);

    std::vector<float> v_match_pt{};
    std::vector<float> v_match_eta{};
    std::vector<float> v_match_phi{};

    v_match_pt.reserve(v_match_pt_size);
    v_match_eta.reserve(v_match_pt_size);
    v_match_phi.reserve(v_match_pt_size);

    auto const& a_match_pt = this->array<float>(matchJetCollBranchName + "_pt");
    auto const& a_match_eta = this->array<float>(matchJetCollBranchName + "_eta");
    auto const& a_match_phi = this->array<float>(matchJetCollBranchName + "_phi");

    for (auto idx = 0; idx < v_match_pt_size; ++idx) {
      float corr = 1;
      if (matchJetCollRequiresJecA) {
        corr = jecA_.correction(a_match_pt[idx], a_match_eta[idx]);
      }

      v_match_pt.emplace_back(a_match_pt[idx] * corr);
      v_match_eta.emplace_back(a_match_eta[idx]);
      v_match_phi.emplace_back(a_match_phi[idx]);
    }

    std::map<size_t, size_t> mapMatchIndices;
    std::vector<float> vecMatchMinDeltaR2(v_pt_size, -1.f);
    for (auto idx : fhDataIndices) {
      int indexBestMatch = -1;
      auto& dR2min = vecMatchMinDeltaR2[idx];
      for (auto idxMatch = 0; idxMatch < v_match_pt_size; ++idxMatch) {
        auto const passesMatchJetPtMin = matchJetPtMin < 0 or v_match_pt[idxMatch] > matchJetPtMin;
        auto const passesMatchJetPtMax = matchJetPtMax < 0 or v_match_pt[idxMatch] < matchJetPtMax;
        if (not(passesMatchJetPtMin and passesMatchJetPtMax)) {
          continue;
        }

        auto const dR2 = utils::deltaR2(v_eta[idx], v_phi[idx], v_match_eta[idxMatch], v_match_phi[idxMatch]);
        if (dR2min < 0 or dR2 < dR2min) {
          dR2min = dR2;
          if (dR2 < matchJetDeltaR2Min) {
            indexBestMatch = idxMatch;
          }
        }
      }

      if (indexBestMatch >= 0) {
        mapMatchIndices.insert(std::make_pair(idx, indexBestMatch));
      }
    }

    // JEC-related histos
    auto const jetCategoryLabelsForJECHistos_v = jetCategoryLabelsForJECHistos(fhData.jetCollection, matchLabel);
    for (auto const& catLabel : jetCategoryLabelsForJECHistos_v) {
      std::vector<size_t> jetIndices;
      jetIndices.reserve(fhDataIndices.size());
      for (auto idx : fhDataIndices) {
        if (jetCategoryForJECFuncMap_[catLabel](std::abs(v_eta[idx]), nCT)) {
          jetIndices.emplace_back(idx);
        }
      }

      for (auto const jetIdx : jetIndices) {
        auto mapMatchIndicesIter(mapMatchIndices.find(jetIdx));
        if (mapMatchIndicesIter == mapMatchIndices.end()) {
          continue;
        }

        auto const jetMatchIdx(mapMatchIndicesIter->first);
        auto const jetMatchPt(v_match_pt[jetMatchIdx]);
        if (jetMatchPt <= 0) {
          continue;
        }

        auto const jetPt(v_pt[jetIdx]);
        auto const jetPtRatio(jetPt / jetMatchPt);
        auto const jetEta(v_eta[jetIdx]);

        H2(dirPrefix + fhData.jetCollection + catLabel + "_MatchedTo" + matchLabel + "_pt_over" + matchLabel + "__vs__" + matchLabel + "_pt")->Fill(jetPtRatio, jetMatchPt, weight);
        H2(dirPrefix + fhData.jetCollection + catLabel + "_MatchedTo" + matchLabel + "_pt_over" + matchLabel + "__vs__pt")->Fill(jetPtRatio, jetPt, weight);
        H1(dirPrefix + fhData.jetCollection + catLabel + "_MatchedTo" + matchLabel + "_eta")->Fill(jetEta, weight);
        H1(dirPrefix + fhData.jetCollection + catLabel + "_MatchedTo" + matchLabel + "_nCT")->Fill(nCT, weight);
        H1(dirPrefix + fhData.jetCollection + catLabel + "_MatchedTo" + matchLabel + "_nPU")->Fill(nPU, weight);
      }
    }
  }
}

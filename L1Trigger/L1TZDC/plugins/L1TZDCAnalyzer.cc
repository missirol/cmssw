//Copied from L1TStage2CaloAnalyzer on 2023.08.04 as it exists in CMSSW_13_1_0_pre4
//Modified by Chris McGinn to instead work for just ZDC etSums
//Contact at christopher.mc.ginn@cern.ch or cfmcginn @ github for bugs

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeter/interface/CaloCluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

//For the output
#include "TTree.h"
//string for some branch handling
#include <string>

//
// class declaration
//

namespace l1t {

  class L1TZDCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit L1TZDCAnalyzer(const edm::ParameterSet&);
    ~L1TZDCAnalyzer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    //fileservice
    edm::Service<TFileService> m_fs;
    //Declare a tree, member and pointer
    TTree* m_zdcEtSumTree_p;
    //Declare the etSum max and bpx max
    //MAX Here is effectively a hardcoded number - you just need to know it from the producer - realistically this is best handled by something passed but for now (2023.08.17) we hardcode
    static const int m_maxBPX = 5;
    float m_zdcEtSumP[m_maxBPX];
    float m_zdcEtSumM[m_maxBPX];

    // ----------member data ---------------------------
    edm::EDGetToken m_sumToken;

    bool doHistos_;

    TFileDirectory evtDispDir_;
  };

  //
  // constants, enums and typedefs
  //

  //
  // static data member definitions
  //

  //
  // constructors and destructor
  //
  L1TZDCAnalyzer::L1TZDCAnalyzer(const edm::ParameterSet& iConfig)
      : doHistos_(iConfig.getUntrackedParameter<bool>("doHistos", true)) {
    usesResource(TFileService::kSharedResource);
    //now do what ever initialization is needed

    // register what you consume and keep token for later access:
    edm::InputTag nullTag("None");

    edm::InputTag sumTag = iConfig.getParameter<edm::InputTag>("etSumToken");
    m_sumToken = consumes<l1t::EtSumBxCollection>(sumTag);

    std::cout << "Processing " << sumTag.label() << std::endl;
  }

  //
  // member functions
  //

  // ------------ method called for each event  ------------
  void L1TZDCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

    //    Handle<EtSumBxCollection> sums;
    Handle<BXVector<l1t::EtSum> > sums;
    iEvent.getByToken(m_sumToken, sums);

    for (int ibx = sums->getFirstBX(); ibx <= sums->getLastBX(); ++ibx) {
      for (auto itr = sums->begin(ibx); itr != sums->end(ibx); ++itr) {
        if (itr->getType() == l1t::EtSum::EtSumType::kZDCP)
          m_zdcEtSumP[ibx] = itr->hwPt();
        if (itr->getType() == l1t::EtSum::EtSumType::kZDCM)
          m_zdcEtSumM[ibx] = itr->hwPt();
      }
    }

    m_zdcEtSumTree_p->Fill();
  }

  // ------------ method called once each job just before starting event loop  ------------
  void L1TZDCAnalyzer::beginJob() {
    m_zdcEtSumTree_p = m_fs->make<TTree>("zdcEtSumTree", "");
    m_zdcEtSumTree_p->Branch("zdcEtSumP", m_zdcEtSumP, ("zdcEtSumP[" + std::to_string(m_maxBPX) + "]/F").c_str());
    m_zdcEtSumTree_p->Branch("zdcEtSumM", m_zdcEtSumM, ("zdcEtSumM[" + std::to_string(m_maxBPX) + "]/F").c_str());
  }

  // ------------ method called once each job just after ending the event loop  ------------
  void L1TZDCAnalyzer::endJob() {}

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void L1TZDCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("etSumToken", edm::InputTag("zdcEtSumProducer", "zdcEtSums"));
    descriptions.add("l1tZDCAnalyzer", desc);
  }

}  // namespace l1t

using namespace l1t;

//define this as a plug-in
DEFINE_FWK_MODULE(L1TZDCAnalyzer);

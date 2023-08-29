// -*- C++ -*-
//
// Package:    L1Trigger/L1TZDC
// Class:      L1TZDC
//
/**\class L1TZDC L1TZDCProducer.cc L1Trigger/L1TZDC/plugins/L1TZDCProducer.cc

 Description: ZDC Producer for L1 Trigger emulation

 Implementation:
     Below text meant to indicate this was largely copied from some original work by James;
modified to be appropriate for producing ZDC Et Sums
*/
//
// Original Author:  James Brooke
//         Created:  Thu, 05 Dec 2013 17:39:27 GMT
//
// Copied for ZDC by: Chris McGinn
//        Copy Made: Wed, 03 Aug 2023
//        Contact: christopher.mc.ginn@cern.ch or
//                 cfmcginn on github for bugs/issues
//
#include <memory>
#include <iostream>
// user include files

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2FirmwareFactory.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2MainProcessor.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsO2ORcd.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

//
// class declaration
//

using namespace l1t;

class L1TZDCProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TZDCProducer(const edm::ParameterSet& ps);
  ~L1TZDCProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //  void endRun(edm::Run const&, edm::EventSetup const&) override;

  int zdcLUTIndexHelper(int iDetPos, int iBXPos);

  // ----------member data ---------------------------

  // input tokens
  //Add the ZDC token, candidateToken for caloParamsHelper
  edm::EDGetTokenT<QIE10DigiCollection> m_zdcToken;
  edm::ESGetToken<CaloParams, L1TCaloParamsRcd> m_candidateToken;

  //input ints
  int bxFirst_;
  int bxLast_;
  int sampleToCenterBX_;

  // put tokens
  edm::EDPutTokenT<EtSumBxCollection> m_etToken;

  //Following the L1TStage2Layer2Producer
  std::unique_ptr<CaloParamsHelper> m_params;
  int m_scaleFactor;
};

L1TZDCProducer::L1TZDCProducer(const edm::ParameterSet& ps) {
  // register what you produce
  m_etToken = produces<EtSumBxCollection>("zdcEtSums");

  // register what you consume and keep token for later access:
  m_zdcToken = consumes<QIE10DigiCollection>(ps.getParameter<edm::InputTag>("zdcToken"));
  m_candidateToken = esConsumes<CaloParams, L1TCaloParamsRcd, edm::Transition::BeginRun>();

  bxFirst_ = ps.getParameter<int>("bxFirst");
  bxLast_ = ps.getParameter<int>("bxLast");
  sampleToCenterBX_ = ps.getParameter<int>("sampleToCenterBX");
}

// ------------ method called to produce the data  ------------
void L1TZDCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace l1t;

  LogDebug("l1t|stage 2") << "L1TZDCProducer::produce function called..." << std::endl;

  // reduced collection to be emplaced in output
  EtSumBxCollection etsumsReduced(0, bxFirst_, bxLast_);

  //inputs
  Handle<QIE10DigiCollection> zdcDigiCollection;
  iEvent.getByToken(m_zdcToken, zdcDigiCollection);

  //Produce ZDC EtSums IF input zdcDigiCollection is valid
  if (zdcDigiCollection.isValid()) {
    //In lieu of bxFirst, bxLast, use the number of the timeslice samples
    const QIE10DataFrame& frametest = (*zdcDigiCollection)[0];
    int nSamples = frametest.samples();

    // Full BX Handling
    int fullBXFirst = -sampleToCenterBX_;
    int fullBXLast = nSamples - sampleToCenterBX_;
    if (nSamples % 2 == 0)
      ++fullBXLast;

    // Full etsum collection, not to be emplaced
    EtSumBxCollection etsums(0, fullBXFirst, fullBXLast);

    //rawadc[detector index][time slices]
    unsigned short rawadc[18][10];

    // the loop below loops over all the elements of the QIE10DigiCollection. Each entry corresponds to one channel
    for (QIE10DigiCollection::const_iterator it = zdcDigiCollection->begin(); it != zdcDigiCollection->end(); it++) {
      const QIE10DataFrame& frame(*it);
      HcalZDCDetId cell = frame.id();
      int zside = cell.zside();
      int section = cell.section();
      int channel = cell.channel();

      if (zside != -1 && zside != 1)
        continue;
      if (section != 1 && section != 2)
        continue;
      if (section == 1 && (channel < 1 || channel > 5))
        continue;
      if (section == 2 && (channel < 1 || channel > 4))
        continue;

      int ihitid = (zside == 1 ? 9 : 0) + (section == 2 ? 5 : 0) + (channel - 1);
      //the loop below iterates over the time slices
      for (int iTS = 0; iTS < nSamples; iTS++) {
        unsigned short adc = (unsigned short)frame[iTS].adc();
        rawadc[ihitid][iTS] = adc;
      }  // end of loop over iTS
    }    //end of loop over channels

    for (int ibx = 0; ibx < nSamples; ibx++) {
      double cEMP = 0, cEMM = 0, cHDP = 0, cHDM = 0;
      double sumcEMP = 0, sumcEMM = 0, sumcHDP = 0, sumcHDM = 0;
      //idet=0-4 correpond to the EM channels
      for (int idet = 0; idet < 5; idet++) {
        unsigned short EMP = rawadc[idet + 9][ibx];
        unsigned short EMM = rawadc[idet][ibx];

        int cEMP_LUTIndex = zdcLUTIndexHelper(idet + 9, (int)EMP);
        int cEMM_LUTIndex = zdcLUTIndexHelper(idet, (int)EMM);

        cEMP = ((double)m_params->zdcLUT()->data(cEMP_LUTIndex)) / ((double)m_scaleFactor);
        cEMM = ((double)m_params->zdcLUT()->data(cEMM_LUTIndex)) / ((double)m_scaleFactor);

        sumcEMP = sumcEMP + cEMP;
        sumcEMM = sumcEMM + cEMM;
      }
      //idet=5-8 correspond to HAD channels
      for (int idet = 5; idet < 9; idet++) {
        unsigned short HDP = rawadc[idet + 9][ibx];
        unsigned short HDM = rawadc[idet][ibx];

        int cHDP_LUTIndex = zdcLUTIndexHelper(idet + 9, (int)HDP);
        int cHDM_LUTIndex = zdcLUTIndexHelper(idet, (int)HDM);

        cHDP = ((double)m_params->zdcLUT()->data(cHDP_LUTIndex)) / ((double)m_scaleFactor);
        cHDM = ((double)m_params->zdcLUT()->data(cHDM_LUTIndex)) / ((double)m_scaleFactor);

        sumcHDP = sumcHDP + cHDP;
        sumcHDM = sumcHDM + cHDM;
      }
      double sumM = sumcEMM + sumcHDM;
      double sumP = sumcEMP + sumcHDP;

      if (ibx == 4) {
        edm::LogInfo("L1TZDCProducer") << ", sumM= " << sumM << std::endl;
        edm::LogInfo("L1TZDCProducer") << ", sumP= " << sumP << std::endl;
      }
      l1t::EtSum tempEtM = l1t::EtSum();
      tempEtM.setHwPt(sumM);
      tempEtM.setHwEta(-1.);
      tempEtM.setHwPhi(0.);
      tempEtM.setType(EtSum::EtSumType::kZDCM);

      l1t::EtSum tempEtP = l1t::EtSum();
      tempEtP.setHwPt(sumP);
      tempEtP.setHwEta(1.);
      tempEtP.setHwPhi(0.);
      tempEtP.setType(EtSum::EtSumType::kZDCP);

      //One object distinguished by type P, M
      etsums.push_back(ibx - sampleToCenterBX_, CaloTools::etSumP4Demux(tempEtP));
      etsums.push_back(ibx - sampleToCenterBX_, CaloTools::etSumP4Demux(tempEtM));

      if (ibx >= sampleToCenterBX_ + bxFirst_ && ibx <= sampleToCenterBX_ + bxLast_) {
        etsumsReduced.push_back(ibx - sampleToCenterBX_, CaloTools::etSumP4Demux(tempEtP));
        etsumsReduced.push_back(ibx - sampleToCenterBX_, CaloTools::etSumP4Demux(tempEtM));
      }
    }  // end of loop over bunch crossings
  }    // end if(zdcDigiCollection.isValid())
  else {
    // If the collection is not valid issue a warning before putting an empty collection
    edm::LogWarning("L1TZDCProducer") << "zdcToken is not valid; return empty ZDC Et Sum BXCollection" << std::endl;
  }

  // Emplace even if !zdcDigiCollection.isValid()
  // Output in this case will be an empty collection
  iEvent.emplace(m_etToken, std::move(etsumsReduced));
}

// ------------ method called when starting to processes a run  ------------
void L1TZDCProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  edm::ESHandle<CaloParams> candidateHandle = iSetup.getHandle(m_candidateToken);

  m_params = std::make_unique<CaloParamsHelper>(*candidateHandle.product());
  m_scaleFactor = m_params->zdcLUT()->data(0);  //First position is the integer scaling factor
  edm::LogInfo("L1TZDCProducer") << "SCALE FACTOR FOR LUT: " << m_scaleFactor << std::endl;
}

// LUT HELPER METHOD
int L1TZDCProducer::zdcLUTIndexHelper(int iDetPos, int iBxPos) { return 1 + iDetPos * 256 + iBxPos; }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TZDCProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("zdcToken", edm::InputTag("hcalDigis", "ZDC"));
  desc.add<int>("bxFirst", -2);
  desc.add<int>("bxLast", 2);
  desc.add<int>("sampleToCenterBX", 4);
  descriptions.add("l1tZDCProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TZDCProducer);

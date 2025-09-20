// -*- C++ -*-
//
// Package:    HLTL1NumberFilter
// Class:      HLTL1NumberFilter
//
/**\class HLTL1NumberFilter HLTL1NumberFilter.cc filter/HLTL1NumberFilter/src/HLTL1NumberFilter.cc

Description:

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Martin Grunewald
//         Created:  Tue Jan 22 13:55:00 CET 2008
//
//

// system include files
#include <string>
#include <iostream>
#include <memory>

// user include files
#include "HLTL1NumberFilter.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/TCDS/interface/TCDSRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// constructors and destructor
//
HLTL1NumberFilter::HLTL1NumberFilter(const edm::ParameterSet& config)
    :  //now do what ever initialization is needed
      inputToken_(consumes<FEDRawDataCollection>(config.getParameter<edm::InputTag>("rawInput"))),
      period_(config.getParameter<unsigned int>("period")),
      fedId_(config.getParameter<int>("fedId")),
      invert_(config.getParameter<bool>("invert")),
      // only try and use TCDS event number if the FED ID 1024 is selected
      useTCDS_(config.getParameter<bool>("useTCDSEventNumber") and fedId_ == 1024) {}

HLTL1NumberFilter::~HLTL1NumberFilter() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

void HLTL1NumberFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("rawInput", edm::InputTag("source"));
  desc.add<unsigned int>("period", 4096);
  desc.add<bool>("invert", true);
  desc.add<int>("fedId", 812);
  desc.add<bool>("useTCDSEventNumber", false);
  descriptions.add("hltL1NumberFilter", desc);
}
//
// member functions
//

// ------------ method called on each new Event  ------------
bool HLTL1NumberFilter::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  using namespace edm;

//  if (iEvent.isRealData()) {
//    bool accept(false);
//    edm::Handle<FEDRawDataCollection> theRaw;
//    iEvent.getByToken(inputToken_, theRaw);
//    const FEDRawData& data = theRaw->FEDData(fedId_);
//    if (data.data() and data.size() > 0) {
//      unsigned long counter;
//      if (useTCDS_) {
//        TCDSRecord record(data.data());
//        counter = record.getTriggerCount();
//      } else {
//        FEDHeader header(data.data());
//        counter = header.lvl1ID();
//      }
//      if (period_ != 0)
//        accept = (counter % period_ == 0);
//      if (invert_)
//        accept = not accept;
//      return accept;
//    } else {
//      LogWarning("HLTL1NumberFilter") << "No valid data for FED " << fedId_ << " used by HLTL1NumberFilter";
//      return false;
//    }
//  } else {
//    return true;
//  }

  auto const evtId{iEvent.id().event()};

  auto const& rawData = iEvent.get(inputToken_);

  auto const& fedData1024 = rawData.FEDData(1024);
  auto const& fedData1230 = rawData.FEDData(1230);
  auto const& fedData1386 = rawData.FEDData(1386);
  auto const& fedData1404 = rawData.FEDData(1404);

  TCDSRecord tcdsRecord{fedData1024.data()};
  auto const tcdsEventNumber{tcdsRecord.getEventNumber()};
  auto const tcdsTriggerCount{tcdsRecord.getTriggerCount()};

  FEDHeader header1024{fedData1024.data()};
  auto const l1Id1024{header1024.lvl1ID()};

  FEDHeader header1230{fedData1230.data()};
  auto const l1Id1230{header1230.lvl1ID()};

  FEDHeader header1386{fedData1386.data()};
  auto const l1Id1386{header1386.lvl1ID()};

  FEDHeader header1404{fedData1404.data()};
  auto const l1Id1404{header1404.lvl1ID()};

//  edm::LogPrint("HLTL1NumberFilter")
//    << " evtId=" << evtId
//    << " tcdsEventNumber=" << tcdsEventNumber
//    << " tcdsTriggerCount=" << tcdsTriggerCount
//    << " l1Id1024=" << l1Id1024
//    << " l1Id1230=" << l1Id1230
//    << " l1Id1386=" << l1Id1386
//    << " l1Id1404=" << l1Id1404
//  ;

  if (evtId != tcdsEventNumber) throw cms::Exception("InputError") << "tcdsEventNumber";

  if (tcdsTriggerCount != l1Id1024) throw cms::Exception("InputError") << "l1Id1024";
  if (tcdsTriggerCount != l1Id1230) throw cms::Exception("InputError") << "l1Id1230";
  if (tcdsTriggerCount != l1Id1386) throw cms::Exception("InputError") << "l1Id1386";
  if (tcdsTriggerCount != l1Id1404) throw cms::Exception("InputError") << "l1Id1404";

  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTL1NumberFilter);

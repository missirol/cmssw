#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtilHelper.h"

l1t::L1TGlobalUtilHelper::L1TGlobalUtilHelper(edm::InputTag const& l1tAlgBlkInputTag,
                                              edm::InputTag const& l1tExtBlkInputTag,
                                              bool const readPrescalesFromFile,
                                              edm::ConsumesCollector& iC)
    : m_l1tAlgBlkInputTag(l1tAlgBlkInputTag),
      m_l1tExtBlkInputTag(l1tExtBlkInputTag),
      m_readPrescalesFromFile(readPrescalesFromFile),
      m_l1tAlgBlkToken(iC.consumes<GlobalAlgBlkBxCollection>(m_l1tAlgBlkInputTag)),
      m_l1tExtBlkToken(iC.consumes<GlobalExtBlkBxCollection>(m_l1tExtBlkInputTag)) {}

l1t::L1TGlobalUtilHelper::L1TGlobalUtilHelper(edm::ParameterSet const& pset, edm::ConsumesCollector& iC)
    : L1TGlobalUtilHelper(pset.getParameter<edm::InputTag>("l1tAlgBlkInputTag"),
                          pset.getParameter<edm::InputTag>("l1tExtBlkInputTag"),
                          pset.getParameter<bool>("ReadPrescalesFromFile"),
                          iC) {}

void l1t::L1TGlobalUtilHelper::fillDescription(edm::ParameterSetDescription& desc,
                                               edm::InputTag const& iAlg,
                                               edm::InputTag const& iExt,
                                               bool const readPrescalesFromFile) {
  desc.add<edm::InputTag>("l1tAlgBlkInputTag", iAlg);
  desc.add<edm::InputTag>("l1tExtBlkInputTag", iExt);
  desc.add<bool>("ReadPrescalesFromFile", readPrescalesFromFile);
}

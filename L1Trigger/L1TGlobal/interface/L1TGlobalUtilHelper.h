#ifndef L1Trigger_L1TGlobal_L1TGlobalUtilHelper_h
#define L1Trigger_L1TGlobal_L1TGlobalUtilHelper_h

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace edm {
  class ConsumesCollector;
  class ParameterSet;
  class ParameterSetDescription;
}

namespace l1t {

  class L1TGlobalUtilHelper {
  public:
    // Using this constructor will require InputTags to be specified in the configuration
    L1TGlobalUtilHelper(edm::ParameterSet const& pset, edm::ConsumesCollector& iC);

    // A module defining its fillDescriptions function might want to use this
    static void fillDescription(edm::ParameterSetDescription& desc,
                                edm::InputTag const& iAlg,
                                edm::InputTag const& iExt,
                                bool readPrescalesFromFile);

    edm::InputTag const& l1tAlgBlkInputTag() const { return m_l1tAlgBlkInputTag; }
    edm::InputTag const& l1tExtBlkInputTag() const { return m_l1tExtBlkInputTag; }

    bool const& readPrescalesFromFile() const { return m_readPrescalesFromFile; }

    edm::EDGetTokenT<GlobalAlgBlkBxCollection> const& l1tAlgBlkToken() const { return m_l1tAlgBlkToken; }
    edm::EDGetTokenT<GlobalExtBlkBxCollection> const& l1tExtBlkToken() const { return m_l1tExtBlkToken; }

  private:
    edm::InputTag const m_l1tAlgBlkInputTag;
    edm::InputTag const m_l1tExtBlkInputTag;

    edm::EDGetTokenT<GlobalAlgBlkBxCollection> const m_l1tAlgBlkToken;
    edm::EDGetTokenT<GlobalExtBlkBxCollection> const m_l1tExtBlkToken;

    bool const m_readPrescalesFromFile;
  };

}  // namespace l1t

#endif

#ifndef L1Trigger_L1TMuonEndCap_PtAssignmentEngineDxy_h
#define L1Trigger_L1TMuonEndCap_PtAssignmentEngineDxy_h

#include <string>

#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "L1Trigger/L1TMuonEndCap/interface/Common.h"
#include "L1Trigger/MLUtilities/interface/HLS4MLModelWrapper.h"

class PtAssignmentEngineAux2017;

class PtAssignmentEngineDxy {
public:
  explicit PtAssignmentEngineDxy(std::string const& modelName);
  virtual ~PtAssignmentEngineDxy() = default;

  void configure(int verbose);

  const PtAssignmentEngineAux2017& aux() const;

  virtual void calculate_pt_dxy(const EMTFTrack& track, emtf::Feature& feature, emtf::Prediction& prediction) const;

  virtual void preprocessing_dxy(const EMTFTrack& track, emtf::Feature& feature) const;

  virtual void call_hls_dxy(const emtf::Feature& feature, emtf::Prediction& prediction) const;

protected:
  l1t::HLS4MLModelWrapper const modelWrapper_;
  int verbose_;
};

#endif

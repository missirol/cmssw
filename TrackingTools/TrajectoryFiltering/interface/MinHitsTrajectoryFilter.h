#ifndef MinHitsTrajectoryFilter_H
#define MinHitsTrajectoryFilter_H

#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

/** A TrajectoryFilter that stops reconstruction if P_t drops
 *  below some value at some confidence level.
 *  The CkfTrajectoryBuilder uses this class to
 *  implement the minimal P_t cut.
 */

class MinHitsTrajectoryFilter final : public TrajectoryFilter {
public:
  explicit MinHitsTrajectoryFilter(int minHits = 5, int seedPairPenalty = 0)
      : theMinHits(minHits), theSeedPairPenalty(seedPairPenalty) {}

  MinHitsTrajectoryFilter(const edm::ParameterSet& pset, edm::ConsumesCollector& iC)
      : theMinHits(pset.getParameter<int>("minimumNumberOfHits")),
        theSeedPairPenalty(pset.getParameter<int>("seedPairPenalty")) {}

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<int>("minimumNumberOfHits", 5);
    iDesc.add<int>("seedPairPenalty", 0);
  }

  bool qualityFilter(const Trajectory& traj) const override { return QF<Trajectory>(traj); }
  bool qualityFilter(const TempTrajectory& traj) const override { return QF<TempTrajectory>(traj); }

  bool toBeContinued(TempTrajectory&) const override { return TrajectoryFilter::toBeContinuedIfNotContributing; }
  bool toBeContinued(Trajectory&) const override { return TrajectoryFilter::toBeContinuedIfNotContributing; }

  std::string name() const override { return "MinHitsTrajectoryFilter"; }

protected:
  template <class T>
  bool QF(const T& traj) const {
    int seedPenalty = (2 == traj.seedNHits()) ? theSeedPairPenalty : 0;  // increase by one if seed-doublet...
    return (traj.foundHits() >= theMinHits + seedPenalty);
  }

  int theMinHits;
  int theSeedPairPenalty;
};

#endif

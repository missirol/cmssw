#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Likely.h"

GeometricSearchDet::~GeometricSearchDet() {}

void GeometricSearchDet::compatibleDetsV(const TrajectoryStateOnSurface& startingState,
                                         const Propagator& prop,
                                         const MeasurementEstimator& est,
                                         std::vector<DetWithState>& result) const {
  if UNLIKELY (!hasGroups())
    edm::LogError("DetLayers") << "At the moment not a real implementation";

  // standard implementation of compatibleDets() for class which have
  // groupedCompatibleDets implemented.

  std::vector<DetGroup> vectorGroups;
edm::LogPrint("GeometricSearchDet") << "GeometricSearchDet::compatibleDetsV-0 " << __LINE__ << " " << startingState.globalPosition() << " " << startingState.localPosition();
  groupedCompatibleDetsV(startingState, prop, est, vectorGroups);
  for (auto itDG = vectorGroups.begin(); itDG != vectorGroups.end(); itDG++) {
    for (auto itDGE = itDG->begin(); itDGE != itDG->end(); itDGE++) {
edm::LogPrint("GeometricSearchDet") << "  GeometricSearchDet::compatibleDetsV-10 " << __LINE__ << " " << itDGE->det() << " " << itDGE->trajectoryState().isValid() << " " << itDGE->trajectoryState().globalPosition();
edm::LogPrint("GeometricSearchDet") << "  GeometricSearchDet::compatibleDetsV-11 " << __LINE__ << " " << itDGE->det()->position() << " " << itDGE->det()->subDetector() << " " << itDGE->det()->geographicalId();
      result.emplace_back(itDGE->det(), itDGE->trajectoryState());
    }
  }
edm::LogPrint("GeometricSearchDet") << "GeometricSearchDet::compatibleDetsV-2 " << __LINE__ << " " << startingState.globalPosition() << " " << startingState.localPosition();
}

void GeometricSearchDet::groupedCompatibleDetsV(const TrajectoryStateOnSurface& startingState,
                                                const Propagator&,
                                                const MeasurementEstimator&,
                                                std::vector<DetGroup>&) const {
  edm::LogError("DetLayers") << "At the moment not a real implementation";
}

std::vector<GeometricSearchDet::DetWithState> GeometricSearchDet::compatibleDets(
    const TrajectoryStateOnSurface& startingState, const Propagator& prop, const MeasurementEstimator& est) const {
  std::vector<DetWithState> result;
edm::LogPrint("GeometricSearchDet") << "GeometricSearchDet::compatibleDets-0 " << __LINE__ << " " << startingState.globalPosition() << " " << startingState.localPosition();
  compatibleDetsV(startingState, prop, est, result);
edm::LogPrint("GeometricSearchDet") << "GeometricSearchDet::compatibleDets-1 " << __LINE__ << " " << startingState.globalPosition() << " " << startingState.localPosition() << " " << result.size();
  return result;
}

std::vector<DetGroup> GeometricSearchDet::groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                                                const Propagator& prop,
                                                                const MeasurementEstimator& est) const {
  std::vector<DetGroup> result;
  groupedCompatibleDetsV(startingState, prop, est, result);
  return result;
}

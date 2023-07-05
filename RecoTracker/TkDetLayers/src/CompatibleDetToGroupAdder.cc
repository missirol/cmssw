#include "CompatibleDetToGroupAdder.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DetGroupMerger.h"

using namespace std;

bool CompatibleDetToGroupAdder::add(const GeometricSearchDet& det,
                                    const TrajectoryStateOnSurface& tsos,
                                    const Propagator& prop,
                                    const MeasurementEstimator& est,
                                    vector<DetGroup>& result) {
  if (det.hasGroups()) {
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    vector<DetGroup> tmp;
    det.groupedCompatibleDetsV(tsos, prop, est, tmp);
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    if (tmp.empty())
      return false;

edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    if (result.empty())
      result.swap(tmp);
    else
      DetGroupMerger::addSameLevel(std::move(tmp), result);
  } else {
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    vector<GeometricSearchDet::DetWithState> compatDets;
    det.compatibleDetsV(tsos, prop, est, compatDets);
    if (compatDets.empty())
      return false;

edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    if (result.empty())
      result.push_back(DetGroup(0, 1));  // empty group for insertion

edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    if (result.size() != 1)
      edm::LogError("TkDetLayers")
          << "CompatibleDetToGroupAdder: det is not grouped but result has more than one group!";
    result.front().reserve(result.front().size() + compatDets.size());
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
    for (vector<GeometricSearchDet::DetWithState>::const_iterator i = compatDets.begin(); i != compatDets.end(); i++)
      result.front().push_back(*i);
  }
  return true;
}

#include "TrackingTools/DetLayers/interface/GeomDetCompatibilityChecker.h"
// #include "TkGeomDetCompatibilityChecker.h"

bool CompatibleDetToGroupAdder::add(const GeomDet& det,
                                    const TrajectoryStateOnSurface& tsos,
                                    const Propagator& prop,
                                    const MeasurementEstimator& est,
                                    vector<DetGroup>& result) {
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
  //TkGeomDetCompatibilityChecker theCompatibilityChecker;
  GeomDetCompatibilityChecker theCompatibilityChecker;
  auto&& compat = theCompatibilityChecker.isCompatible(&det, tsos, prop, est);
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
  if (!compat.first)
    return false;
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
  if (result.empty())
    result.push_back(DetGroup(0, 1));  // empty group for ge insertion

edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
  if (result.size() != 1)
    edm::LogError("TkDetLayers") << "CompatibleDetToGroupAdder: det is not grouped but result has more than one group!";

edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
  result.front().emplace_back(&det, std::move(compat.second));
edm::LogPrint("CompatibleDetToGroupAdder") << " CompatibleDetToGroupAdder::add " << __LINE__ << " tsos.globalPos=" << tsos.globalPosition();
  return true;
}

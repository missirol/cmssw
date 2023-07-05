#ifndef RecHitPropagator_H
#define RecHitPropagator_H

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class TrackingRecHit;
class MagneticField;
class Plane;

class dso_hidden RecHitPropagator {
public:
  TrajectoryStateOnSurface propagate(const TrackingRecHit& hit,
                                     const Plane& plane,
                                     const TrajectoryStateOnSurface& ts) const;
};

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

// propagate from glued to mono/stereo
inline TrajectoryStateOnSurface fastProp(const TrajectoryStateOnSurface& ts, const Plane& oPlane, const Plane& tPlane) {
  GlobalVector gdir = ts.globalMomentum();
edm::LogPrint("RecHitPropagator") << " RecHitPropagator::fastProp " << __LINE__ << " " << gdir << " " << oPlane.position();
  double delta = tPlane.localZ(oPlane.position());
edm::LogPrint("RecHitPropagator") << " RecHitPropagator::fastProp " << __LINE__ << " " << gdir << " " << tPlane.position() << " " << delta;
  LocalVector ldir = tPlane.toLocal(gdir);  // fast prop!
edm::LogPrint("RecHitPropagator") << " RecHitPropagator::fastProp " << __LINE__ << " " << ldir << " " << ts.globalPosition();
  LocalPoint lPos = tPlane.toLocal(ts.globalPosition());
  LocalPoint projectedPos = lPos - ldir * delta / ldir.z();
edm::LogPrint("RecHitPropagator") << " RecHitPropagator::fastProp " << __LINE__ << " " << lPos << " " << projectedPos;
  // we can also patch it up as only the position-errors are used...
  GlobalTrajectoryParameters gp(
      tPlane.toGlobal(projectedPos), gdir, ts.charge(), &ts.globalParameters().magneticField());
  if (ts.hasError())
    return TrajectoryStateOnSurface(gp, ts.curvilinearError(), tPlane);
  else
    return TrajectoryStateOnSurface(gp, tPlane);
}

#endif

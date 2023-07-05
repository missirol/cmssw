#include "TrackingTools/GeomPropagators/interface/StraightLinePropagator.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "TrackingTools/GeomPropagators/interface/PropagationExceptions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

std::pair<TrajectoryStateOnSurface, double> StraightLinePropagator::propagateWithPath(const FreeTrajectoryState& fts,
                                                                                      const Plane& plane) const {
  // propagate parameters
  LocalPoint x;
  LocalVector p;
  double s = 0;
  bool parametersOK = propagateParametersOnPlane(fts, plane, x, p, s);
edm::LogPrint("StraightLinePropagator") << "   StraightLinePropagator::propagateWithPath " << __LINE__ << " parametersOK=" << parametersOK
<< " x=" << x
<< " p=" << p
<< " s=" << s
;
  if (!parametersOK)
    return std::make_pair(TrajectoryStateOnSurface(), 0.);

  // compute propagated state
  if (fts.hasError()) {
    return std::make_pair(propagatedState(fts, plane, jacobian(s), x, p), s);
  } else {
    // return state without errors
    return std::make_pair(TSOS(LocalTrajectoryParameters(x, p, fts.charge()), plane, theField), s);
  }
}

std::pair<TrajectoryStateOnSurface, double> StraightLinePropagator::propagateWithPath(const FreeTrajectoryState& fts,
                                                                                      const Cylinder& cylinder) const {
  // propagate parameters
  GlobalPoint x;
  GlobalVector p;
  double s = 0;
  bool parametersOK = propagateParametersOnCylinder(fts, cylinder, x, p, s);
  if (!parametersOK)
    return std::make_pair(TrajectoryStateOnSurface(), 0.);

  // compute propagated state
  if (fts.hasError()) {
    return std::make_pair(propagatedState(fts, cylinder, jacobian(s), x, p), s);
  } else {
    // return state without errors
    return std::make_pair(TSOS(GlobalTrajectoryParameters(x, p, fts.charge(), theField), cylinder), s);
  }
}

TrajectoryStateOnSurface StraightLinePropagator::propagatedState(const FTS& fts,
                                                                 const Surface& surface,
                                                                 const AlgebraicMatrix55& jacobian,
                                                                 const LocalPoint& x,
                                                                 const LocalVector& p) const {
  if (fts.hasError()) {
    // propagate error
    TSOS tmp(fts, surface);
    const AlgebraicSymMatrix55& eLocal = tmp.localError().matrix();
    AlgebraicSymMatrix55 lte = ROOT::Math::Similarity(jacobian, eLocal);
    LocalTrajectoryError eloc(lte);
    LocalTrajectoryParameters ltp(x, p, fts.charge());
    return TSOS(ltp, eloc, surface, theField);
  } else {
    // return state without errors
    return TSOS(LocalTrajectoryParameters(x, p, fts.charge()), surface, theField);
  }
}

TrajectoryStateOnSurface StraightLinePropagator::propagatedState(const FTS& fts,
                                                                 const Surface& surface,
                                                                 const AlgebraicMatrix55& jacobian,
                                                                 const GlobalPoint& x,
                                                                 const GlobalVector& p) const {
  if (fts.hasError()) {
    // propagate error
    TSOS tmp(fts, surface);
    const AlgebraicSymMatrix55& eLocal = tmp.localError().matrix();
    AlgebraicSymMatrix55 lte = ROOT::Math::Similarity(jacobian, eLocal);
    LocalTrajectoryError eloc(lte);

    TSOS tmp2(tmp.localParameters(), eloc, surface, theField);
    GlobalTrajectoryParameters gtp(x, p, fts.charge(), theField);
    return TSOS(gtp, tmp2.cartesianError(), surface);
  } else {
    // return state without errors
    return TSOS(GlobalTrajectoryParameters(x, p, fts.charge(), theField), surface);
  }
}

AlgebraicMatrix55 StraightLinePropagator::jacobian(double& s) const {
  //Jacobian for 5*5 local error matrix
  AlgebraicMatrix55 j = AlgebraicMatrixID();  //Jacobian

  double dir = (propagationDirection() == alongMomentum) ? 1. : -1.;
  if (s * dir < 0.)
    return j;

  j(3, 1) = s;
  j(4, 2) = s;

  return j;
}

bool StraightLinePropagator::propagateParametersOnCylinder(
    const FTS& fts, const Cylinder& cylinder, GlobalPoint& x, GlobalVector& p, double& s) const {
  GlobalPoint sp = cylinder.toGlobal(LocalPoint(0., 0.));
  if (sp.x() != 0. || sp.y() != 0.) {
    throw PropagationException("Cannot propagate to an arbitrary cylinder");
  }

  x = fts.position();
  p = fts.momentum();
  s = cylinder.radius() - x.perp();

  double dir = (propagationDirection() == alongMomentum) ? 1. : -1.;
  if (s * dir < 0.)
    return false;

  double dt = s / p.perp();
  x = GlobalPoint(x.x() + p.x() * dt, x.y() + p.y() * dt, x.z() + p.z() * dt);

  return true;
}

bool StraightLinePropagator::propagateParametersOnPlane(
    const FTS& fts, const Plane& plane, LocalPoint& x, LocalVector& p, double& s) const {
edm::LogPrint("StraightLinePropagator") << "   StraightLinePropagator::propagateParametersOnPlane " << __LINE__ << " x=" << x << " p=" << p << " s=" << s << " fts.pos=" << fts.position() << " fts.mom=" << fts.momentum();
  //Do extrapolation in local frame of plane
  //  LocalPoint sp = plane.toLocal(plane.position());


//   StraightLinePropagator::propagateParametersOnPlane 130 x= (0,0,0)  p= (0,0,0)  s=0 fts.pos= (0.116712,-0.182284,-1.7058)  fts.mom= (33.2922,19.2505,23.3234)
//   StraightLinePropagator::propagateParametersOnPlane 142 x= (-22.7081,-28.586,12.9057)  p= (30.8843,32.697,-9.53674e-07)  s=-12.9057
//   StraightLinePropagator::propagateWithPath 17 parametersOK=1 x= (4.17944e+08,4.42474e+08,0)  p= (30.8843,32.697,-9.53674e-07)  s=-12.9057
//    GeomDetCompatibilityChecker::isCompatible 123 geoId=369141860 tsos.globalPos= (0.116712,-0.182284,-1.7058)  propSt.globalPos= (4.50529e+08,2.60508e+08,3.15626e+08)
//    GeomDetCompatibilityChecker::isCompatible 134 geoId=369141860 tsos.globalPos= (0.116712,-0.182284,-1.7058)  propSt.globalPos= (4.50529e+08,2.60508e+08,3.15626e+08)  es=1


  x = plane.toLocal(fts.position());
  p = plane.toLocal(fts.momentum());
  s = -x.z();  // sp.z() - x.z(); local z of plane always 0

  //double dir = (propagationDirection() == alongMomentum) ? 1. : -1.;
  //if(s*dir < 0.) return false;
  if ((p.x() != 0 || p.y() != 0) && p.z() == 0 && s != 0)
    return false;

edm::LogPrint("StraightLinePropagator") << "   StraightLinePropagator::propagateParametersOnPlane " << __LINE__ << " x=" << x << " p=" << p << " s=" << s;

  x = LocalPoint(x.x() + (p.x() / p.z()) * s, x.y() + (p.y() / p.z()) * s, x.z() + s);

  return true;
}

/*
 * \class TrackerSeedCleaner
 *  Reference class for seeds cleaning
 *  Seeds Cleaner based on sharedHits cleaning, direction cleaning and pt cleaning
    \author A. Grelli -  Purdue University, Pavia University
 */

#include "RecoMuon/TrackerSeedGenerator/interface/TrackerSeedCleaner.h"

//---------------
// C++ Headers --
//---------------
#include <vector>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include "RecoTracker/TkTrackingRegions/interface/TkTrackingRegionsMargin.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoRange.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

using namespace std;
using namespace edm;

//
// inizialization
//
void TrackerSeedCleaner::init(const MuonServiceProxy* service) {
  theProxyService = service;
}

//
//
//
void TrackerSeedCleaner::setEvent(const edm::Event& event) { event.getByToken(beamspotToken_, bsHandle_); }

//
// clean seeds
//
void TrackerSeedCleaner::clean(const reco::TrackRef& muR,
                               const RectangularEtaPhiTrackingRegion& region,
                               tkSeeds& seeds) {
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
  // call the shared input cleaner
  if (cleanBySharedHits)
    seeds = define(seeds);
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
  theTTRHBuilder = theProxyService->eventSetup().getHandle(theTTRHBuilderToken);
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

  LogDebug("TrackerSeedCleaner") << seeds.size() << " trajectory seeds to the events before cleaning" << endl;

  //check the validity otherwise vertexing
  const reco::BeamSpot& bs = *bsHandle_;
  /*reco track and seeds as arguments. Seeds eta and phi are checked and 
   based on deviation from L2 eta and phi seed is accepted or not*/

edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
  std::vector<TrajectorySeed> result;

edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
  TSCBLBuilderNoMaterial tscblBuilder;
  // PerigeeConversions tspConverter;
  for (TrajectorySeedCollection::iterator seed = seeds.begin(); seed < seeds.end(); ++seed) {
    if (seed->nHits() < 2)
      continue;
    //get parameters and errors from the seed state
    TransientTrackingRecHit::RecHitPointer recHit = theTTRHBuilder->build(&*(seed->recHits().end() - 1));
    TrajectoryStateOnSurface state = trajectoryStateTransform::transientState(
        seed->startingState(), recHit->surface(), theProxyService->magneticField().product());
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    TrajectoryStateClosestToBeamLine tsAtClosestApproachSeed =
        tscblBuilder(*state.freeState(), bs);  //as in TrackProducerAlgorithms
    if (!tsAtClosestApproachSeed.isValid())
      continue;
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
    GlobalPoint vSeed1 = tsAtClosestApproachSeed.trackStateAtPCA().position();
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
    GlobalVector pSeed = tsAtClosestApproachSeed.trackStateAtPCA().momentum();
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
    GlobalPoint vSeed(vSeed1.x() - bs.x0(), vSeed1.y() - bs.y0(), vSeed1.z() - bs.z0());
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    //eta,phi info from seeds
    double etaSeed = state.globalMomentum().eta();
    double phiSeed = pSeed.phi();

    //if the limits are too stringent rescale limits
    typedef PixelRecoRange<float> Range;
    typedef TkTrackingRegionsMargin<float> Margin;

    Range etaRange = region.etaRange();
    double etaLimit = (fabs(fabs(etaRange.max()) - fabs(etaRange.mean())) < 0.1)
                          ? 0.1
                          : fabs(fabs(etaRange.max()) - fabs(etaRange.mean()));
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    Margin phiMargin = region.phiMargin();
    double phiLimit = (phiMargin.right() < 0.1) ? 0.1 : phiMargin.right();

    double ptSeed = pSeed.perp();
    double ptMin = (region.ptMin() > 3.5) ? 3.5 : region.ptMin();
    // Clean
    bool inEtaRange = etaSeed >= (etaRange.mean() - etaLimit) && etaSeed <= (etaRange.mean() + etaLimit);
    bool inPhiRange = (fabs(deltaPhi(phiSeed, double(region.direction().phi()))) < phiLimit);
    // pt cleaner
    bool inPtRange = ptSeed >= ptMin && ptSeed <= 2 * (muR->pt());
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    // save efficiency don't clean triplets with pt cleaner
    if (seed->nHits() == 3)
      inPtRange = true;
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    // use pt and angle cleaners
    if (inPtRange && usePt_Cleaner && !useDirection_Cleaner) {
      result.push_back(*seed);
      LogDebug("TrackerSeedCleaner") << " Keeping the seed : this seed passed pt selection";
    }
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    // use only angle default option
    if (inEtaRange && inPhiRange && !usePt_Cleaner && useDirection_Cleaner) {
      result.push_back(*seed);
      LogDebug("TrackerSeedCleaner") << " Keeping the seed : this seed passed direction selection";
    }
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    // use all the cleaners
    if (inEtaRange && inPhiRange && inPtRange && usePt_Cleaner && useDirection_Cleaner) {
      result.push_back(*seed);
      LogDebug("TrackerSeedCleaner") << " Keeping the seed : this seed passed direction and pt selection";
    }
edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;

    LogDebug("TrackerSeedCleaner") << " eta for current seed " << etaSeed << "\n"
                                   << " phi for current seed " << phiSeed << "\n"
                                   << " eta for L2 track  " << muR->eta() << "\n"
                                   << " phi for L2 track  " << muR->phi() << "\n";
  }

edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
  //the new seeds collection
  if (!result.empty() && (useDirection_Cleaner || usePt_Cleaner))
    seeds.swap(result);

edm::LogPrint("AAA") << "TrackerSeedCleaner " << __LINE__;
  LogDebug("TrackerSeedCleaner") << seeds.size() << " trajectory seeds to the events after cleaning" << endl;

  return;
}





//
// the sharedHits cleaner
//

std::vector<TrajectorySeed> TrackerSeedCleaner::define(std::vector<TrajectorySeed> const& coll) const {

  std::vector<TrajectorySeed> result;
  result.reserve(coll.size());

  std::vector<bool> maskPairs(coll.size(), true);

  for (size_t i1 = 0; i1 < coll.size(); ++i1) {
    auto const& s1 = coll[i1];
    if (s1.nHits() == 3) {
      continue;
    }

    for (size_t i2 = i1 + 1; i2 < coll.size(); ++i2) {
      auto const& s2 = coll[i2];
      if (s2.nHits() != 3) {
        continue;
      }

      int shared = 0;
      for (auto const& h2 : s2.recHits()) {
        for (auto const& h1 : s1.recHits()) {
          if (h2.sharesInput(&h1, TrackingRecHit::all))
            ++shared;
        }
      }

      if (shared == 2) {
        maskPairs[i1] = false;
        break;
      }
    }
  }

  for (size_t i1 = 0; i1 < coll.size(); ++i1) {
    if (maskPairs[i1]) {
      result.emplace_back(coll[i1]);
    }
  }

  return result;
}

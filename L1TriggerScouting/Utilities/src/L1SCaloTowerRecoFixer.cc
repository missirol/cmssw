#include <cassert>
#include <cmath>

#include "FWCore/Utilities/interface/Exception.h"

#include "L1TriggerScouting/Utilities/interface/L1SCaloTowerRecoFixer.h"

L1SCaloTowerRecoFixer::L1SCaloTowerRecoFixer(edm::RunNumber_t const run, edm::EventNumber_t const orbit, int const bx) {
  assert(orbit > 0);
  assert(bx > 0);

  orbit_correct_ = orbit_input_ = orbit;
  bx_correct_ = bx_input_ = bx;

  if (run < 401681) {
    orbit_correct_ -= 1;
    bx_correct_ -= 2;
  } else if (401681 <= run and run < 403259) {
    orbit_correct_ -= 1;
    bx_correct_ -= 3;
  } else if (403462 <= run and run < 403552) {
    orbit_correct_ -= 1;
    bx_correct_ -= 2;
  } else if (403655 <= run and run < 403681) {
    bx_correct_ -= 2;
  }

  if (bx_correct_ < 1) {
    bx_correct_ += 3564;

    assert(orbit_correct_ > 0);
    orbit_correct_ -= 1;
  }

  assert(bx_correct_ > 0 and bx_correct_ < 3565);
  assert(orbit_correct_ > 0);

  isFromMP70_ = ((bx_correct_ % 9) == 7);
  isFromMP71_ = ((bx_correct_ % 9) == 8);

  // Runs prior to the online deployment of https://github.com/cms-sw/cmssw/pull/50469
  shiftPhiByPlus1_ = (run < 402330);

  // https://gitlab.cern.ch/scouting-demonstrator/calol2/-/merge_requests/13
  applyEtaPhiFix1_ = (run < 403259 or (403462 <= run and run < 403552));

  // https://gitlab.cern.ch/scouting-demonstrator/calol2/-/merge_requests/14
  applyEtaPhiFix2_ = (run < 403655) and isFromMP71_;
}

L1SCaloTowerRecoFixer::CaloTowerHwEtaAndHwPhi L1SCaloTowerRecoFixer::correctCaloTowerHwEtaAndHwPhi(
    int16_t const hwEta, int16_t const hwPhi) const {
  int16_t const old_hwEta = hwEta;
  int16_t const old_hwEtaAbs = std::abs(old_hwEta);
  int16_t const old_hwPhi = shiftPhiByPlus1_ ? (hwPhi + 1) : hwPhi;

  if (old_hwEtaAbs < 1 or old_hwEtaAbs > 28) {
    throw cms::Exception("InvalidValue") << "invalid value of old_hwEta: " << old_hwEta;
  }

  if (old_hwPhi < 1 or old_hwPhi > 72) {
    throw cms::Exception("InvalidValue") << "invalid value of old_hwPhi: " << old_hwPhi;
  }

  if (not(applyEtaPhiFix1_ or applyEtaPhiFix2_)) {
    return {old_hwEta, old_hwPhi};
  }

  bool isFirstWord = false;
  int16_t old_link = -1;

  if (applyEtaPhiFix1_) {
    isFirstWord = (old_hwEta > 0);
    old_link = (old_hwPhi - 1);
  } else {
    isFirstWord = ((old_hwPhi % 2) != 0);
    if (old_hwEta < 0) {
      old_link = ((old_hwPhi % 2) == 0) ? (old_hwPhi - 1) : old_hwPhi;
    } else {
      old_link = ((old_hwPhi % 2) == 0) ? (old_hwPhi - 2) : (old_hwPhi - 1);
    }
  }

  int16_t new_link = old_link;
  if (applyEtaPhiFix2_) {
    auto newLinkIt = kMP71LinkMap.find(old_link);
    if (newLinkIt != kMP71LinkMap.end()) {
      new_link = newLinkIt->second;
    }
  }

  int16_t const new_hwEta = ((new_link % 2) != 0) ? -old_hwEtaAbs : old_hwEtaAbs;
  int16_t const new_hwEtaAbs = std::abs(new_hwEta);

  int16_t const link_phi = ((new_link % 2) != 0) ? (new_link - 1) : new_link;
  int16_t const new_hwPhi = isFirstWord ? (link_phi + 1) : (link_phi + 2);

  if (new_hwEtaAbs < 1 or new_hwEtaAbs > 28) {
    throw cms::Exception("InvalidValue") << "invalid value of new_hwEta: " << new_hwEta;
  }

  if (new_hwPhi < 1 or new_hwPhi > 72) {
    throw cms::Exception("InvalidValue") << "invalid value of new_hwPhi: " << new_hwPhi;
  }

  return {new_hwEta, new_hwPhi};
}

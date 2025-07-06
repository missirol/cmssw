#ifndef DataFormats_Scouting_Run3ScoutingEERecHit_h
#define DataFormats_Scouting_Run3ScoutingEERecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for EERecHits
//
// IMPORTANT: any changes to Run3ScoutingEERecHit must be backward-compatible !

class Run3ScoutingEERecHit {
public:
  Run3ScoutingEERecHit(float energy, int ix, int iy, bool positiveZ)
      : energy_{energy}, ix_{ix}, iy_{iy}, positiveZ_{positiveZ} {}

  Run3ScoutingEERecHit() : energy_{0}, ix_{0}, iy_{0}, positiveZ_{false} {}

  float energy() const { return energy_; }
  int ix() const { return ix_; }
  int iy() const { return iy_; }
  bool positiveZ() { return positiveZ_; }

private:
  float energy_;
  int ix_;
  int iy_;
  bool positiveZ_;
};

using Run3ScoutingEERecHitCollection = std::vector<Run3ScoutingEERecHit>;

#endif

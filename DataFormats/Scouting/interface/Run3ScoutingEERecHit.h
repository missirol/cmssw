#ifndef DataFormats_Scouting_Run3ScoutingEERecHit_h
#define DataFormats_Scouting_Run3ScoutingEERecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for EERecHits
//
// IMPORTANT: any changes to Run3ScoutingEERecHit must be backward-compatible !

class Run3ScoutingEERecHit {
public:
  Run3ScoutingEERecHit(float energy, uint8_t ix, uint8_t iy, bool positiveZ)
      : energy_{energy}, ix_{ix}, iy_{iy}, positiveZ_{positiveZ} {}

  Run3ScoutingEERecHit() : energy_{0}, ix_{0}, iy_{0}, positiveZ_{false} {}

  float energy() const { return energy_; }
  uint8_t ix() const { return ix_; }
  uint8_t iy() const { return iy_; }
  bool positiveZ() { return positiveZ_; }

private:
  float energy_;
  uint8_t ix_;
  uint8_t iy_;
  bool positiveZ_;
};

using Run3ScoutingEERecHitCollection = std::vector<Run3ScoutingEERecHit>;

#endif

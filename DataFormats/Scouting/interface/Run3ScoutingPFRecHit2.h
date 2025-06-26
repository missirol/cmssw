#ifndef DataFormats_Scouting_Run3ScoutingPFRecHit2_h
#define DataFormats_Scouting_Run3ScoutingPFRecHit2_h

#include <vector>

// Run-3 HLT-Scouting data format for PFRecHits
//
// IMPORTANT: any changes to Run3ScoutingPFRecHit2 must be backward-compatible !

class Run3ScoutingPFRecHit2 {
public:
  Run3ScoutingPFRecHit2(uint32_t detId, float energy) : detId_{detId}, energy_{energy} {}

  Run3ScoutingPFRecHit2() : detId_{0}, energy_{0} {}

  unsigned int detId() const { return detId_; }
  float energy() const { return energy_; }

private:
  unsigned int detId_;
  float energy_;
};

using Run3ScoutingPFRecHit2Collection = std::vector<Run3ScoutingPFRecHit2>;

#endif

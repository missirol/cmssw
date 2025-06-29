#ifndef DataFormats_Scouting_Run3ScoutingPFRecHit_h
#define DataFormats_Scouting_Run3ScoutingPFRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for PFRecHits
//
// IMPORTANT: any changes to Run3ScoutingPFRecHit must be backward-compatible !

class Run3ScoutingPFRecHit {
public:
  Run3ScoutingPFRecHit(unsigned int detId, float energy) : detId_{detId}, energy_{energy} {}

  Run3ScoutingPFRecHit() : detId_{0}, energy_{0} {}

  unsigned int detId() const { return detId_; }
  float energy() const { return energy_; }

private:
  unsigned int detId_;
  float energy_;
};

using Run3ScoutingPFRecHitCollection = std::vector<Run3ScoutingPFRecHit>;

#endif

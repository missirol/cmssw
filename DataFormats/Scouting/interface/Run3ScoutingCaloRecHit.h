#ifndef DataFormats_Scouting_Run3ScoutingCaloRecHit_h
#define DataFormats_Scouting_Run3ScoutingCaloRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for CaloRecHits
//
// IMPORTANT: any changes to Run3ScoutingCaloRecHit must be backward-compatible !

class Run3ScoutingCaloRecHit {
public:
  Run3ScoutingCaloRecHit(float energy, unsigned int detId) : energy_{energy}, detId_{detId} {}

  Run3ScoutingCaloRecHit() : energy_{0}, detId_{0} {}

  float energy() const { return energy_; }
  unsigned int detId() const { return detId_; }

private:
  float energy_;
  unsigned int detId_;
};

using Run3ScoutingCaloRecHitCollection = std::vector<Run3ScoutingCaloRecHit>;

#endif

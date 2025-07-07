#ifndef DataFormats_Scouting_Run3ScoutingCaloRecHit_h
#define DataFormats_Scouting_Run3ScoutingCaloRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for CaloRecHits
//
// IMPORTANT: any changes to Run3ScoutingCaloRecHit must be backward-compatible !

class Run3ScoutingCaloRecHit {
public:
  Run3ScoutingCaloRecHit(float energy, unsigned int detId, uint32_t flags) : energy_{energy}, detId_{detId}, flags_{flags} {}

  Run3ScoutingCaloRecHit() : energy_{0}, detId_{0}, flags_{0} {}

  float energy() const { return energy_; }
  unsigned int detId() const { return detId_; }
  uint32_t flags() const { return flags_; }

private:
  float energy_;
  unsigned int detId_;
  uint32_t flags_;
};

using Run3ScoutingCaloRecHitCollection = std::vector<Run3ScoutingCaloRecHit>;

#endif

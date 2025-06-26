#ifndef DataFormats_Scouting_Run3ScoutingCaloRecHit_h
#define DataFormats_Scouting_Run3ScoutingCaloRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for calorimeter RecHits
//
// IMPORTANT: any changes to Run3ScoutingCaloRecHit must be backward-compatible !

class Run3ScoutingCaloRecHit {
public:
  Run3ScoutingCaloRecHit(uint32_t id, float energy) : id_{id}, energy_{energy} {}

  Run3ScoutingCaloRecHit() : id_{0}, energy_{0} {}

  uint32_t id() const { return id_; }
  float energy() const { return energy_; }

private:
  uint32_t id_;
  float energy_;
};

using Run3ScoutingCaloRecHitCollection = std::vector<Run3ScoutingCaloRecHit>;

#endif

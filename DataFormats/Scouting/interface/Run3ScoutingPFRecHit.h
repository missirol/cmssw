#ifndef DataFormats_Scouting_Run3ScoutingPFRecHit_h
#define DataFormats_Scouting_Run3ScoutingPFRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for PFRecHits
//
// IMPORTANT: any changes to Run3ScoutingPFRecHit must be backward-compatible !

class Run3ScoutingPFRecHit {
public:
  Run3ScoutingPFRecHit(uint32_t detId, float energy, float rho, float eta, float phi) : detId_{detId}, energy_{energy}, rho_{rho}, eta_{eta}, phi_{phi} {}

  Run3ScoutingPFRecHit() : detId_{0}, energy_{0}, rho_{0}, eta_{0}, phi_{0} {}

  unsigned int detId() const { return detId_; }
  float energy() const { return energy_; }
  float rho() const { return rho_; }
  float eta() const { return eta_; }
  float phi() const { return phi_; }

private:
  unsigned int detId_;
  float energy_;
  float rho_;
  float eta_;
  float phi_;
};

using Run3ScoutingPFRecHitCollection = std::vector<Run3ScoutingPFRecHit>;

#endif

#ifndef DataFormats_Scouting_Run3ScoutingPFRecHit_h
#define DataFormats_Scouting_Run3ScoutingPFRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for PFRecHits
//
// IMPORTANT: any changes to Run3ScoutingPFRecHit must be backward-compatible !

class Run3ScoutingPFRecHit {
public:
  Run3ScoutingPFRecHit(float energy, float rho, float eta, float phi) : energy_{energy}, rho_{rho}, eta_{eta}, phi_{phi} {}

  Run3ScoutingPFRecHit() : energy_{0}, rho_{0}, eta_{0}, phi_{0} {}

  float energy() const { return energy_; }
  float rho() const { return rho_; }
  float eta() const { return eta_; }
  float phi() const { return phi_; }

private:
  float energy_;
  float rho_;
  float eta_;
  float phi_;
};

using Run3ScoutingPFRecHitCollection = std::vector<Run3ScoutingPFRecHit>;

#endif

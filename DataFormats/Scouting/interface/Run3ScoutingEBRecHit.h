#ifndef DataFormats_Scouting_Run3ScoutingEBRecHit_h
#define DataFormats_Scouting_Run3ScoutingEBRecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for EBRecHits
//
// IMPORTANT: any changes to Run3ScoutingEBRecHit must be backward-compatible !

class Run3ScoutingEBRecHit {
public:
  Run3ScoutingEBRecHit(float energy, int ieta, int iphi) : energy_{energy}, ieta_{ieta}, iphi_{iphi} {}

  Run3ScoutingEBRecHit() : energy_{0}, ieta_{0}, iphi_{0} {}

  float energy() const { return energy_; }
  int ieta() const { return ieta_; }
  int iphi() const { return iphi_; }

private:
  float energy_;
  int ieta_;
  int iphi_;
};

using Run3ScoutingEBRecHitCollection = std::vector<Run3ScoutingEBRecHit>;

#endif

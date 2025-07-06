#ifndef DataFormats_Scouting_Run3ScoutingHBHERecHit_h
#define DataFormats_Scouting_Run3ScoutingHBHERecHit_h

#include <vector>

// Run-3 HLT-Scouting data format for HBHERecHits
//
// IMPORTANT: any changes to Run3ScoutingHBHERecHit must be backward-compatible !

class Run3ScoutingHBHERecHit {
public:
  Run3ScoutingHBHERecHit(float energy, int ieta, int iphi, int depth)
      : energy_{energy}, ieta_{ieta}, iphi_{iphi}, depth_{depth} {}

  Run3ScoutingHBHERecHit() : energy_{0}, ieta_{0}, iphi_{0}, depth_{0} {}

  float energy() const { return energy_; }
  int ieta() const { return ieta_; }
  int iphi() const { return iphi_; }
  int depth() const { return depth_; }

private:
  float energy_;
  int ieta_;
  int iphi_;
  int depth_;
};

using Run3ScoutingHBHERecHitCollection = std::vector<Run3ScoutingHBHERecHit>;

#endif

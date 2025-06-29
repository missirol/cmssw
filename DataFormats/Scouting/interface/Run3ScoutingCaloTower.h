#ifndef DataFormats_Scouting_Run3ScoutingCaloTower_h
#define DataFormats_Scouting_Run3ScoutingCaloTower_h

#include <vector>

// Run-3 HLT-Scouting data format corresponding to
// DataFormats/CaloTowers/interface/CaloTower.h
//
// IMPORTANT: any changes to Run3ScoutingCaloTower must be backward-compatible !

class Run3ScoutingCaloTower {
public:
  Run3ScoutingCaloTower(float emEnergy, float hadEnergy, float eta, float phi)
      : emEnergy_{emEnergy}, hadEnergy_{hadEnergy}, eta_{eta}, phi_{phi} {}

  Run3ScoutingCaloTower() : emEnergy_{0}, hadEnergy_{0}, eta_{0}, phi_{0} {}

  float emEnergy() const { return emEnergy_; }
  float hadEnergy() const { return hadEnergy_; }
  float eta() const { return eta_; }
  float phi() const { return phi_; }

private:
  float emEnergy_;
  float hadEnergy_;
  float eta_;
  float phi_;
};

using Run3ScoutingCaloTowerCollection = std::vector<Run3ScoutingCaloTower>;

#endif

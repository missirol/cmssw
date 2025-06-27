#ifndef DataFormats_Scouting_Run3ScoutingCaloTower2_h
#define DataFormats_Scouting_Run3ScoutingCaloTower2_h

#include <vector>

// Run-3 HLT-Scouting data format corresponding to
// DataFormats/CaloTowers/interface/CaloTower.h
//
// IMPORTANT: any changes to Run3ScoutingCaloTower2 must be backward-compatible !

class Run3ScoutingCaloTower2 {
public:
  Run3ScoutingCaloTower2(float emEnergy,
                        float hadEnergy,
                        float outerEnergy,
                        float eta,
                        float phi)
      : emEnergy_{emEnergy},
        hadEnergy_{hadEnergy},
        outerEnergy_{outerEnergy},
        eta_{ieta},
        phi_{iphi} {}

  Run3ScoutingCaloTower2()
      : emEnergy_{0},
        hadEnergy_{0},
        outerEnergy_{0},
        eta_{0},
        phi_{0} {}

  float emEnergy() const { return emEnergy_; }
  float hadEnergy() const { return hadEnergy_; }
  float outerEnergy() const { return outerEnergy_; }
  float eta() const { return eta_; }
  float phi() const { return phi_; }

private:
  float emEnergy_;
  float hadEnergy_;
  float outerEnergy_;
  float eta_;
  float phi_;
};

using Run3ScoutingCaloTower2Collection = std::vector<Run3ScoutingCaloTower2>;

#endif

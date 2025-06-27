#ifndef DataFormats_Scouting_Run3ScoutingCaloTower_h
#define DataFormats_Scouting_Run3ScoutingCaloTower_h

#include <vector>

// Run-3 HLT-Scouting data format corresponding to
// DataFormats/CaloTowers/interface/CaloTower.h
//
// IMPORTANT: any changes to Run3ScoutingCaloTower must be backward-compatible !

class Run3ScoutingCaloTower {
public:
  Run3ScoutingCaloTower(float pt,
                        float eta,
                        float phi,
                        float m,
                        float emEnergy,
                        float hadEnergy,
                        float outerEnergy,
                        float ecalTime,
                        float hcalTime,
                        int ieta,
                        int iphi,
                        int numCrystals,
                        unsigned int numConstituents,
                        uint32_t towerStatusWord)
      : pt_{pt},
        eta_{eta},
        phi_{phi},
        m_{m},
        emEnergy_{emEnergy},
        hadEnergy_{hadEnergy},
        outerEnergy_{outerEnergy},
        ecalTime_{ecalTime},
        hcalTime_{hcalTime},
        ieta_{ieta},
        iphi_{iphi},
        numCrystals_{numCrystals},
        numConstituents_{numConstituents},
        towerStatusWord_{towerStatusWord} {}

  Run3ScoutingCaloTower()
      : pt_{0},
        eta_{0},
        phi_{0},
        m_{0},
        emEnergy_{0},
        hadEnergy_{0},
        outerEnergy_{0},
        ecalTime_{0},
        hcalTime_{0},
        ieta_{0},
        iphi_{0},
        numCrystals_{0},
        numConstituents_{0},
        towerStatusWord_{0} {}

  float pt() const { return pt_; }
  float eta() const { return eta_; }
  float phi() const { return phi_; }
  float m() const { return m_; }

  float emEnergy() const { return emEnergy_; }
  float hadEnergy() const { return hadEnergy_; }
  float outerEnergy() const { return outerEnergy_; }

  float ecalTime() const { return ecalTime_; }
  float hcalTime() const { return hcalTime_; }

  int ieta() const { return ieta_; }
  int iphi() const { return iphi_; }

  int numCrystals() const { return numCrystals_; }
  unsigned int numConstituents() const { return numConstituents_; }

  uint32_t towerStatusWord() const { return towerStatusWord_; }

private:
  float pt_;
  float eta_;
  float phi_;
  float m_;

  float emEnergy_;
  float hadEnergy_;
  float outerEnergy_;

  float ecalTime_;
  float hcalTime_;

  int ieta_;
  int iphi_;

  int numCrystals_;
  unsigned int numConstituents_;

  uint32_t towerStatusWord_;
};

using Run3ScoutingCaloTowerCollection = std::vector<Run3ScoutingCaloTower>;

#endif

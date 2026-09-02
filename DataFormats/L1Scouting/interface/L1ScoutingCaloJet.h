#ifndef DataFormats_L1Scouting_L1ScoutingCaloJet_h
#define DataFormats_L1Scouting_L1ScoutingCaloJet_h

#include <cstdint>

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"

namespace l1ScoutingRun3 {

  class CaloJet {
  public:
    CaloJet()
        : pt_(0), eta_(0), phi_(0), mass_(0), energyCorr_(1), energyFracEm_(0), nConst_(0), nConstSaturatedEnergy_(0) {}

    CaloJet(float pt,
            float eta,
            float phi,
            float mass,
            float energyCorr,
            float energyFracEm,
            int nConst,
            uint32_t nConstSaturatedEnergy)
        : pt_(pt),
          eta_(eta),
          phi_(phi),
          mass_(mass),
          energyCorr_(energyCorr),
          energyFracEm_(energyFracEm),
          nConst_(nConst),
          nConstSaturatedEnergy_(nConstSaturatedEnergy) {}

    CaloJet(float pt,
            float eta,
            float phi,
            float mass,
            float energyCorr,
            float energyFracEm,
            int nConst,
            uint32_t nConstSatEnergyECAL,
            uint32_t nConstSatEnergyHCAL,
            uint32_t nConstSatEnergyBoth)
        : pt_(pt),
          eta_(eta),
          phi_(phi),
          mass_(mass),
          energyCorr_(energyCorr),
          energyFracEm_(energyFracEm),
          nConst_(nConst) {
      setNConstSaturatedEnergy(nConstSatEnergyECAL, nConstSatEnergyHCAL, nConstSatEnergyBoth);
    }

    void setPt(float pt) { pt_ = pt; }
    void setEta(float eta) { eta_ = eta; }
    void setPhi(float phi) { phi_ = phi; }
    void setMass(float mass) { mass_ = mass; }
    void setEnergyCorr(float energyCorr) { energyCorr_ = energyCorr; }
    void setEnergyFracEm(float energyFracEm) { energyFracEm_ = energyFracEm; }
    void setNConst(int nConst) { nConst_ = nConst; }
    void setNConstSaturatedEnergy(uint32_t nConstSaturatedEnergy) { nConstSaturatedEnergy_ = nConstSaturatedEnergy; }
    void setNConstSaturatedEnergy(uint32_t nECAL, uint32_t nHCAL, uint32_t nBoth) {
      nConstSaturatedEnergy_ = 0;
      nConstSaturatedEnergy_ |= (nECAL & 0x7FF);
      nConstSaturatedEnergy_ |= (nHCAL & 0x7FF) << 11;
      nConstSaturatedEnergy_ |= (nBoth & 0x3FF) << 22;
    }

    float pt() const { return pt_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float mass() const { return mass_; }
    float energyCorr() const { return energyCorr_; }
    float energyFracEm() const { return energyFracEm_; }
    int nConst() const { return nConst_; }
    uint32_t nConstSaturatedEnergy() const { return nConstSaturatedEnergy_; }
    uint32_t nConstSaturatedEnergyECAL() const { return nConstSaturatedEnergy_ & 0x7FF; }
    uint32_t nConstSaturatedEnergyHCAL() const { return (nConstSaturatedEnergy_ >> 11) & 0x7FF; }
    uint32_t nConstSaturatedEnergyECALAndHCAL() const { return (nConstSaturatedEnergy_ >> 22) & 0x3FF; }
    uint32_t nConstSaturatedEnergyECALOrHCAL() const {
      return (nConstSaturatedEnergyECAL() + nConstSaturatedEnergyHCAL() - nConstSaturatedEnergyECALAndHCAL());
    }

  private:
    // Jet 4-momentum (pT,eta,phi,mass) after jet energy-scale correction
    float pt_;
    float eta_;
    float phi_;
    float mass_;

    // Jet energy-scale correction factor applied to the jet's 4-momentum
    float energyCorr_;

    // Electromagnetic fraction of the jet's total energy (as measured in ECAL)
    float energyFracEm_;

    // Number of jet constituents (constituent == CaloTower)
    int nConst_;

    // Value encoding 3 non-negative numbers.
    //  - Bits  0-10 (11 bits): number of jet constituents with saturated energy in ECAL.
    //  - Bits 11-21 (11 bits): number of jet constituents with saturated energy in HCAL.
    //  - Bits 22-31 (10 bits): number of jet constituents with saturated energy in both ECAL and HCAL.
    // The number of jet constituents with saturated energy in either ECAL or HCAL
    // can be derived from the 3 numbers above as "ECAL + HCAL - Both".
    uint32_t nConstSaturatedEnergy_;
  };

  using CaloJetOrbitCollection = OrbitCollection<CaloJet>;

}  // namespace l1ScoutingRun3
#endif  // DataFormats_L1Scouting_L1ScoutingCaloJet_h

#ifndef L1TriggerScouting_Utilities_L1SCaloTowerRecoFixer_h
#define L1TriggerScouting_Utilities_L1SCaloTowerRecoFixer_h

#include <cstdint>
#include <unordered_map>

#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"

class L1SCaloTowerRecoFixer {
public:
  explicit L1SCaloTowerRecoFixer(edm::RunNumber_t const run, edm::EventNumber_t const orbit, int const bx);

  edm::EventNumber_t orbit_input() const { return orbit_input_; }
  edm::EventNumber_t orbit_correct() const { return orbit_correct_; }

  int bx_input() const { return bx_input_; }
  int bx_correct() const { return bx_correct_; }

  bool isFromMP70() const { return isFromMP70_; }
  bool isFromMP71() const { return isFromMP71_; }

  bool shiftPhiByPlus1() const { return shiftPhiByPlus1_; }
  bool applyEtaPhiFix1() const { return applyEtaPhiFix1_; }
  bool applyEtaPhiFix2() const { return applyEtaPhiFix2_; }

  struct CaloTowerHwEtaAndHwPhi {
    int16_t hwEta = 0;
    int16_t hwPhi = 0;
  };

  CaloTowerHwEtaAndHwPhi correctCaloTowerHwEtaAndHwPhi(int16_t const hwEta, int16_t const hwPhi) const;

private:
  inline static const std::unordered_map<int16_t, int16_t> kMP71LinkMap = {
      {24, 35},
      {25, 32},
      {26, 34},
      {27, 33},
      {28, 30},
      {29, 31},
      {30, 28},
      {31, 29},
      {32, 25},
      {33, 27},
      {34, 26},
      {35, 24},
  };

  edm::EventNumber_t orbit_input_;
  edm::EventNumber_t orbit_correct_;

  int bx_input_;
  int bx_correct_;

  bool isFromMP70_;
  bool isFromMP71_;

  bool shiftPhiByPlus1_;
  bool applyEtaPhiFix1_;
  bool applyEtaPhiFix2_;
};

#endif

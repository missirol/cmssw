#include <iostream>

#include "L1Trigger/MLUtilities/interface/HLS4MLModelWrapper.h"

int main(int, char**) {
  l1t::HLS4MLModelWrapper mw_a{"GTADModel_v1"};
  l1t::HLS4MLModelWrapper mw_b{"GTADModel_v1"};
  l1t::HLS4MLModelWrapper mw_c{"GTADModel_v3"};

  mw_a.reset();
  mw_b.reset("GTADModel_v3");

  mw_a = std::move(mw_c);
  mw_a.reset("GTADModel_v3");
  mw_a.reset("GTADModel_v1");

  return 0;
}

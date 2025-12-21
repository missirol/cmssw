#include <iostream>

#include "L1Trigger/MLUtilities/interface/HLS4MLModelWrapper.h"

int main(int, char**) {
  l1t::HLS4MLModelWrapper mw1{"GTADModel_v1"};
  l1t::HLS4MLModelWrapper mw3{"GTADModel_v3"};
  l1t::HLS4MLModelWrapper mw4{"GTADModel_v4"};
  l1t::HLS4MLModelWrapper mw5{"GTADModel_v5"};

  mw1.reset();
  mw1.reset("GTADModel_v1");
  mw1.reset("GTADModel_v5");

  return 0;
}

#include "ap_fixed.h"
#include <bitset>
#include <iostream>
#include <string>

template<class T>
std::bitset<18> print(T arg) {
  uint32_t u;
  std::memcpy(&u, &arg, sizeof(u));
  return std::bitset<18>(u);
}

int main(int argc, char** argv) {

  int const a = std::stoi(argv[1]);
  double const b = a * 0.5;

  ap_fixed<18, 13> ret1 = b;

  ap_ufixed<12,11> ret2_tmp = b;
  ap_fixed<18, 13> ret2 = ret2_tmp;

  std::cout << "ret1 = " << ret1 << std::endl;
  std::cout << "ret2 = " << ret2 << std::endl;

  auto const b1 = print(ret1);
  auto const b2 = print(ret2);

  std::cout << "b1 = " << b1 << std::endl;
  std::cout << "b2 = " << b2 << std::endl;
  std::cout << (b1 == b2) << std::endl;

  std::cout << "---------------" << std::endl;

  return 0;
}

#ifndef CommonTools_PileupAlgos_PuppiCandidates_h
#define CommonTools_PileupAlgos_PuppiCandidates_h

#include <array>

template<unsigned int N>
class PuppiCandidates {
 public:
  unsigned int size() const { return N; }
  std::array<float,N> pt;
  std::array<float,N> eta;
  std::array<float,N> rapidity;
  std::array<float,N> phi;
  std::array<int,N> id;
};

#endif

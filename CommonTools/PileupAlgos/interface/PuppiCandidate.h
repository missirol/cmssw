#ifndef CommonTools_PileupAlgos_PuppiCandidate_h
#define CommonTools_PileupAlgos_PuppiCandidate_h

#include <cstdint>

struct PuppiCandidate {
 public:
  PuppiCandidate(const float Pt=0.f, const float Eta=0.f, const float Phi=0.f, const int8_t Id=-1)
    : pt(Pt), eta(Eta), phi(Phi), id(Id) {}
  float pt;
  float eta;
  float phi;
  int8_t id;
};

#endif

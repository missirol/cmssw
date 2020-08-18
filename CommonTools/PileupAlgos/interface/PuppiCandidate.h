#ifndef CommonTools_PileupAlgos_PuppiCandidate_h
#define CommonTools_PileupAlgos_PuppiCandidate_h

class PuppiCandidate {
public:
  PuppiCandidate() : pt(0), eta(0), phi(0), m(0), rapidity(0), id(0) {}
  ~PuppiCandidate() {}

  double pt, eta, phi, m, rapidity;
  int id;
};
#endif

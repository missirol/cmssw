#ifndef CommonTools_PileupAlgos_PuppiCandidate_h
#define CommonTools_PileupAlgos_PuppiCandidate_h

struct PuppiCandidate {
public:
  PuppiCandidate() {
    pt = 0.f;
    eta = 0.f;
    phi = 0.f;
    mass = 0.f;
    time = 0.f;
    depth = 0.f;
    dZ = 0.f;
    d0 = 0.f;
    charge = 0;
    pdgId = 0;
    id = -1;
  }

  float pt;
  float eta;
  float phi;
  float mass;
  float time;
  float depth;
  float dZ;
  float d0;
  int charge;
  int pdgId;
  int id;
};

#endif

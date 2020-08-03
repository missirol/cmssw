#ifndef CommonTools_PileupAlgos_PuppiCandidate_h
#define CommonTools_PileupAlgos_PuppiCandidate_h

struct PuppiCandidate {
public:
  PuppiCandidate() {
    pt = 0.;
    eta = 0.;
    phi = 0.;
    mass = 0.;
    time = 0.;
    depth = 0.;
    dZ = 0.;
    d0 = 0.;
    charge = 0;
    pdgId = 0;
    id = -1;
  }

  double pt;
  double eta;
  double phi;
  double mass;
  double time;
  double depth;
  double dZ;
  double d0;
  int charge;
  int pdgId;
  int id;
};

#endif

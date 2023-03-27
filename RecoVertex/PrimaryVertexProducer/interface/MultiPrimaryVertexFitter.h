#ifndef MultiPrimaryVertexFitter_h
#define MultiPrimaryVertexFitter_h

/**\class MultiPrimaryVertexFitter
 
  Description: simultaneaous fit of primary vertices

*/

//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"

typedef ROOT::Math::SMatrix<double,3> Error3;

class MultiPrimaryVertexFitter : public PrimaryVertexFitterBase {
 public:
  //MultiPrimaryVertexFitter(const edm::ParameterSet &conf);
  MultiPrimaryVertexFitter(double chi2cutoff=2.5, double mintrkweight=0.2);
   ~MultiPrimaryVertexFitter() override = default;

  std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack> &, const std::vector<TransientVertex> &, const reco::BeamSpot &, const bool) override;

 protected:
  std::vector<reco::TransientTrack> input_tracks;
  
  struct TrackInfo{
    float ipsig;// temp
    float x,y;
    float z;
    float odz2;
    double S11, S22, S12;
    Error3 C;
    double c[3];
    double a1[3],a2[3];
    double b1,b2;
    double d;
    //unsigned int kmin, kmax;
    std::vector<double> weight;
  };
  
  void fill_trackinfo(const std::vector<reco::TransientTrack> &, const reco::BeamSpot &);
  void clean(const std::vector<reco::TransientTrack> & tracks, const std::vector<TransientVertex> & clusters);
  void fill_weights(const double beta, const reco::BeamSpot &, const double Zcutoff=0.);
  void dump(const std::string &, const::std::vector<TransientVertex> & clusters, const std::vector<reco::TransientTrack> & tracks, const reco::BeamSpot & beamspot, const double zmin, const double zmax, const unsigned int nit);
  double single_fit(const reco::BeamSpot &, float beam_weight, const bool fill_covariances = false);
  std::vector<TrackInfo> trackinfo;

  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  std::vector<double> rho_vtx;
  std::vector<Error3> V_vtx;
  std::vector<float> chi2_vtx;
  std::vector<float> d2Fperp_vtx;

  int hdump=0; // FIXME debugging only

  // configuration
  double chi2_cutoff_;
  double min_trackweight_;

};
#endif

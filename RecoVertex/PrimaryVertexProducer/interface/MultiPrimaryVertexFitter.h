#ifndef RecoVertex_PrimaryVertexProducer_MultiPrimaryVertexFitter_h
#define RecoVertex_PrimaryVertexProducer_MultiPrimaryVertexFitter_h

/**\class MultiPrimaryVertexFitter

  Description: simultaneaous fit of primary vertices

*/
#include <vector>

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"

class MultiPrimaryVertexFitter : public PrimaryVertexFitterBase {
public:
  //MultiPrimaryVertexFitter(const edm::ParameterSet &conf);
  MultiPrimaryVertexFitter(double const chi2cutoff = 2.5, double const mintrkweight = 0.2, bool const verbose = false);
  ~MultiPrimaryVertexFitter() override = default;

  std::vector<TransientVertex> fit(std::vector<reco::TransientTrack> const&,
                                   std::vector<TransientVertex> const&,
                                   reco::BeamSpot const&,
                                   bool const) override;

protected:
  using Error3 = ROOT::Math::SMatrix<double, 3>;

  std::vector<reco::TransientTrack> input_tracks;

  struct TrackInfo {
    float ipsig;  // temp
    float x, y;
    float z;
    float odz2;
    double S11, S22, S12;
    Error3 C;
    double c[3];
    double a1[3], a2[3];
    double b1, b2;
    double d;
    //unsigned int kmin, kmax;
    std::vector<double> weight;
  };

  void fill_trackinfo(const std::vector<reco::TransientTrack> &, const reco::BeamSpot &);
  void clean(const std::vector<reco::TransientTrack> &tracks, const std::vector<TransientVertex> &clusters);
  void fill_weights(const double beta, const reco::BeamSpot &, const double Zcutoff = 0.);
  void dump(const std::string &,
            const ::std::vector<TransientVertex> &clusters,
            const std::vector<reco::TransientTrack> &tracks,
            const reco::BeamSpot &beamspot,
            const double zmin,
            const double zmax,
            const unsigned int nit);
  double single_fit(const reco::BeamSpot &, float beam_weight, const bool fill_covariances = false);
  std::vector<TrackInfo> trackinfo;

  void test_chisquared(int const k, double const xb, double const yb, double const zb, TrackInfo const& ti) const;

  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  std::vector<double> rho_vtx;
  std::vector<Error3> V_vtx;
  std::vector<float> chi2_vtx;
  std::vector<float> d2Fperp_vtx;

  int hdump = 0;  // FIXME debugging only

  // configuration
  double chi2_cutoff_;
  double min_trackweight_;
  const bool verbose_ = false;
};

#endif

#ifndef RecoVertex_PrimaryVertexProducer_AMultiPrimaryVertexFitter_h
#define RecoVertex_PrimaryVertexProducer_AMultiPrimaryVertexFitter_h

/**\class MultiPrimaryVertexFitter

  Description: simultaneous chisquared fit of primary vertices 

*/
#include <vector>

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"


class AMultiPrimaryVertexFitter : public PrimaryVertexFitterBase {
public:
  AMultiPrimaryVertexFitter(double chicutoff = 2.5, double zcutoff = 1.0, double mintrkweight = 0.4);
  ~AMultiPrimaryVertexFitter() override = default;

  std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack> &,
                                   const std::vector<TransientVertex> &,
                                   const reco::BeamSpot &,
                                   const bool) override;

  using Error3 = ROOT::Math::SMatrix<double, 3>;

protected:

  std::vector<reco::TransientTrack> input_tracks;

  struct TrackInfo {
    double S11, S22, S12;     // inverse of the covariance (sub-)matrix
    Error3 C;                 // HT S H
    double g[3];
    double H1[3], H2[3];
    double b1, b2;
    double zpca, dzError;
  };

  double track_in_vertex_chsq(const TrackInfo &, const double, const double, const double);
  void fill_trackinfo(const std::vector<reco::TransientTrack> &, const reco::BeamSpot &);
  void fill_weights(const double beta, const reco::BeamSpot &);
  Error3 get_inverse_beam_covariance(const reco::BeamSpot &);
  void make_vtx_trk_map(const double beta);
  double update(const reco::BeamSpot &, float beam_weight, const bool fill_covariances = false);
  void remove_vertex(unsigned int);

  void verify(){
    unsigned int nt = trackinfo.size();
    unsigned int nv = xv.size();
    assert( (yv.size() == nv) && "yv size");
    assert( (zv.size() == nv) && "zv size");
    assert( (rho_vtx.size() == nv) && "rho_vtx size");
    assert( (tkfirstv.size() == (nv+1)) && "tkfirstv size");
    assert( (tkmap.size() == tkweight.size()) && "tkmapsize <> tkweightssize");

    for(unsigned int k = 0; k< nv; k++){
      assert( (tkfirstv[k] < tkweight.size()) && "tkfirst[k]");
      assert( (tkfirstv[k+1] < tkweight.size()) && "tkfirst[k+1]");
      assert( (tkfirstv[k] <= tkfirstv[k+1]) && "tkfirst[k/k+1]");
      for(unsigned int j = tkfirstv[k]; j < tkfirstv[k+1]; k++){
	assert( (j<tkmap.size()) && "illegal tkfirst entry");
	unsigned int i = tkmap[j];
	assert( (i<nt) && "illegal tkmap entry");
	assert( (tkweight[i] >=0 ) && "negative tkweight or nan");
	assert( (tkweight[i] <=1 ) && "tkweight > 1 or nan");
      }
    }

    double rho_sum = 0;
    for(unsigned int k = 0; k< nv; k++){
      assert( (rho_vtx[k]==rho_vtx[k]) && "nan in rho_vtx");
      assert( (rho_vtx[k] >= 0) && "negative rho or nan");
      assert( (rho_vtx[k] <= 1) && "rho > 1 or nan");
      rho_sum += rho_vtx[k];
    }
    //assert( (rho_sum <= 1) && "rho_sum > 1");

  };

  // track information
  std::vector<TrackInfo> trackinfo;

  // vertex lists
  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> zv;
  std::vector<Error3> V_vtx;
  std::vector<double> rho_vtx;
  // track-vertex-mapping and weights after a coarse z-cut
  std::vector<unsigned int> tkfirstv;  // parallel to the vertex list
  std::vector<unsigned int> tkmap;     // parallel to tkweight
  std::vector<double> tkweight;        // parallel to tkmap
  std::vector<double> tkw_c; // parallel to the track list, may not be needed
  // configuration
  double chi_cutoff_;
  double z_cutoff_;
  double min_trackweight_;
};
#endif

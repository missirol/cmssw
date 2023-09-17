#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/AMultiPrimaryVertexFitter.h"
#include "vdt/vdtMath.h"

//#define PVTX_DEBUG

AMultiPrimaryVertexFitter::AMultiPrimaryVertexFitter(const double chicutoff,
                                                     const double zcutoff,
                                                     const double mintrkweight)
    : chi_cutoff_(chicutoff), z_cutoff_(zcutoff), min_trackweight_(mintrkweight) {
  edm::LogInfo("AMultiPrimaryVertexFitter") << "instantiating MPV with "
                                            << " chi cutoff =" << chi_cutoff_ << " z cutoff =" << z_cutoff_
                                            << " min_trackweight=" << min_trackweight_ << std::endl;
}

double AMultiPrimaryVertexFitter::track_in_vertex_chsq(const TrackInfo &ti,
                                                       const double xvtx,
                                                       const double yvtx,
                                                       const double zvtx) {
  double F1 = ti.b1 + xvtx * ti.H1[0] + yvtx * ti.H1[1];         // using H1[2]=0
  double F2 = ti.b2 + xvtx * ti.H2[0] + yvtx * ti.H2[1] - zvtx;  // using H2[2]=-1
  double chsq = F1 * F1 * ti.S11 + F2 * F2 * ti.S22 + 2. * F1 * F2 * ti.S12;
#ifdef PVTX_DEBUG
  assert((chsq >= 0) && " negative chi**2");
#endif
  return chsq;
}

void AMultiPrimaryVertexFitter::fill_trackinfo(const std::vector<reco::TransientTrack> &tracks,
                                               const reco::BeamSpot &beamSpot) {
  /* fill track information used during fits into arrays, parallell to the list of input tracks */

  trackinfo.clear();
  trackinfo.reserve(tracks.size());

  for (auto &trk : tracks) {
    TrackInfo ti;
    // F1,F2 can be the perigee parameters (3,4)
    const auto tspca = trk.stateAtBeamLine().trackStateAtPCA();  // freeTrajectoryState
    const auto tspca_pe = PerigeeConversions::ftsToPerigeeError(tspca);
    const auto momentum = tspca.momentum();
    auto const cos_phi = momentum.x() / momentum.perp();
    auto const sin_phi = momentum.y() / momentum.perp();
    auto const tan_lambda = momentum.z() / momentum.perp();

    // covariance matrix of (F1,F2)
    double V11 = tspca_pe.covarianceMatrix()(3, 3);
    double V22 = tspca_pe.covarianceMatrix()(4, 4);
    double V12 = tspca_pe.covarianceMatrix()(3, 4);

    // S = V^{-1}
    double DetV = V11 * V22 - V12 * V12;
    if (fabs(DetV) < 1.e-16) {
      edm::LogWarning("AMultiPrimaryVertexFitter")
          << "Warning, det(V) almost vanishes : " << DetV << " !! This should not happen!" << std::endl;
      ti.S11 = 0;
      ti.S22 = 0;
      ti.S12 = 0;
    } else {
      ti.S11 = V22 / DetV;
      ti.S22 = V11 / DetV;
      ti.S12 = -V12 / DetV;
    }
#ifdef PVTX_DEBUG
    assert((ti.S11 >= 0) && "S11 positivity test");
    assert((ti.S22 >= 0) && "S22 positivity test");
#endif
    // chi**2(track;vertex) = F^T V F
    // F1 = -(x-x0) * sin(phi) + (y-y0)* cos(phi)
    // F2 = z0  + (x-x0)*cos(phi) + (y-y0)*sin(phi) * tan(lambda)- z
    //
    // F_{1,2} (x) = H_{1,2} dot [x,y,z] + b_{1,2}
    // F1 =-x * sin_phi             +y * cos_phi            + 0   + (  x0*sin_phi - y0*cos_phi )
    // F2 = x * cos_phi*tan_lambda + y * sin_phi*tan_lambda - z   + (-(x0*cos_phi + y0*sin_phi)*tan_lambda + z0)
    //
    // the signs of F1,F2 matter for the correlation term
    ti.b1 = tspca.position().x() * sin_phi - tspca.position().y() * cos_phi;
    ti.H1[0] = -sin_phi;
    ti.H1[1] = cos_phi;
    ti.H1[2] = 0;
    ti.b2 = tspca.position().z() - (tspca.position().x() * cos_phi + tspca.position().y() * sin_phi) * tan_lambda;
    ti.H2[0] = cos_phi * tan_lambda;
    ti.H2[1] = sin_phi * tan_lambda;
    ti.H2[2] = -1.;

    // g = b^T S H  [2^T x (2,2) x (2,3) -> 3]     = 2^T x (2,3)
    // C = H^T S H  [(3,2) x (2,2) x (2,3) -> 3x3] = (3,2) x ( 2,3) (symmetric)
    for (int k = 0; k < 3; k++) {
      double SH1k = (ti.S11 * ti.H1[k] + ti.S12 * ti.H2[k]);
      double SH2k = (ti.S12 * ti.H1[k] + ti.S22 * ti.H2[k]);
      ti.g[k] = ti.b1 * SH1k + ti.b2 * SH2k;
      for (int l = 0; l < 3; l++) {
        ti.C(l, k) = ti.H1[l] * SH1k + ti.H2[l] * SH2k;
      }
    }

    ti.zpca = tspca.position().z();
    ti.dzError = trk.track().dzError();
    trackinfo.push_back(ti);
  }
}

void AMultiPrimaryVertexFitter::make_vtx_trk_map(double beta) {
  unsigned const int nv = xv.size();
  unsigned const int nt = trackinfo.size();

  if (nv < 1) {
    edm::LogWarning("AMultiPrimaryVertexFit") << " empty vertex list with " << nt << " tracks" << std::endl;
    return;
  }

  // parallel lists for track to vertex mapping
  // each vertex has a pointer into these lists (tkfirstv)
  tkfirstv.clear();
  tkmap.clear();
  tkweight.clear();

  if (nv == 1) {
    // always accept all tracks for a single vertex fit
    tkfirstv.push_back(0);
    tkfirstv.push_back(nt);
    tkweight.assign(nt, 0.);
    tkmap.reserve(nt);
    for (unsigned int i = 0; i < nt; i++) {
      tkmap.emplace_back(i);
    }
    return;
  }

  // n > 1
  tkmap.reserve(nv * 100);
  tkweight.reserve(nv * 100);
  for (unsigned int k = 0; k < nv; k++) {
    tkfirstv.emplace_back(tkmap.size());
    for (unsigned int i = 0; i < nt; i++) {
      auto &ti = trackinfo[i];
      const double zrange = 5 * sqrt(beta) * ti.dzError;
      if (std::abs(zv[k] - ti.zpca) < z_cutoff_) {
        const double dztrk = ti.b2 + xv[k] * ti.H2[0] + yv[k] * ti.H2[1] - zv[k];
        if (std::abs(dztrk) < zrange) {
          tkmap.emplace_back(i);
          tkweight.emplace_back(1.);
        }
      }
    }
  }
  tkfirstv.push_back(tkmap.size());  // extra entry, simplifies loops, every vertex has a "successor" now
}

void AMultiPrimaryVertexFitter::fill_weights(double beta, const reco::BeamSpot &beamspot) {
  unsigned const int nt = trackinfo.size();
  unsigned const int nv = xv.size();

  double Z_cutoff = vdt::fast_exp(-0.5 * beta * chi_cutoff_ * chi_cutoff_) / nv;

  std::vector<double> Z_track(nt, Z_cutoff);
  const double argmax = std::max(16., 0.5 * beta * chi_cutoff_ * chi_cutoff_ * 5);

  // evaluate and cache track-vertex assignment chi**2 for all clusters and sum up Z
  for (unsigned int k = 0; k < nv; k++) {
    for (unsigned int j = tkfirstv[k]; j < tkfirstv[k + 1]; j++) {
      const unsigned int i = tkmap[j];
      double arg = 0.5 * beta * track_in_vertex_chsq(trackinfo[i], xv[k], yv[k], zv[k]);
      if (arg < argmax) {
        const double e = rho_vtx[k] * vdt::fast_exp(-arg);
        tkweight[j] = e;  // normalized later by the proper Z_track[i]
        Z_track[i] += e;  // sum up exponentials to normalize
      } else {
        tkweight[j] = 0.;
      }
    }
  }

  // now we have the partition function, Z_i and can evaluate assignment probabilities (weights)
  std::vector<double> sumtk_w(nt, 0);
  for (unsigned int j = 0; j < tkmap.size(); j++) {
    const unsigned int i = tkmap[j];
    if (Z_track[i] > 0) {
      tkweight[j] /= Z_track[i];
      sumtk_w[i] += tkweight[j];
    } else {
      tkweight[j] = 0;
    }
  }

  // update vertex masses
  for (unsigned int k = 0; k < nv; k++) {
    rho_vtx[k] = 0;
    for (unsigned int j = tkfirstv[k]; j < tkfirstv[k + 1]; j++) {
      if (sumtk_w[tkmap[j]] > 0) {
        rho_vtx[k] += tkweight[j] / sumtk_w[tkmap[j]];
      }
#ifdef PVTX_DEBUG
      else {
        std::cout << "track weight vanishing for all vertices  nr" << tkmap[j] "/" >> nt << std::endl;
      }
#endif
    }
    rho_vtx[k] /= nt;
  }
}

void AMultiPrimaryVertexFitter::remove_vertex(unsigned int k) {
  unsigned int nv = xv.size();
  if (nv < 2)
    return;

  // give the rho_vtx to the nearest neighbour
  if (k == 0) {
    rho_vtx[1] += rho_vtx[0];
  } else if (k == nv) {
    rho_vtx[k - 1] += rho_vtx[k];
  } else if (std::abs(zv[k + 1] - zv[k]) < std::abs(zv[k] = zv[k - 1])) {
    rho_vtx[k + 1] += rho_vtx[k];
  } else {
    rho_vtx[k - 1] += rho_vtx[k];
  }

  // 1) remove the vertex from the vertex list
  xv.erase(xv.begin() + k);
  yv.erase(yv.begin() + k);
  zv.erase(zv.begin() + k);
  rho_vtx.erase(rho_vtx.begin() + k);
  V_vtx.erase(V_vtx.begin() + k);

  // 2) adjust the track-map map
  // 2a) remove the map entries that belong the deleted vertex
  const unsigned int num_erased_map_entries = tkfirstv[k + 1] - tkfirstv[k];
  tkmap.erase(tkmap.begin() + tkfirstv[k], tkmap.begin() + tkfirstv[k + 1]);
  tkweight.erase(tkweight.begin() + tkfirstv[k], tkweight.begin() + tkfirstv[k + 1]);
  // 2b) adjust pointers for the following vertices, including the dummy entry behind the last (now [nv-1])
  for (unsigned int k1 = k + 1; k1 < nv + 1; k1++) {
    tkfirstv[k1] -= num_erased_map_entries;
  }
  // 2c) erase the pointer of the removed vertex
  tkfirstv.erase(tkfirstv.begin() + k);
}

double AMultiPrimaryVertexFitter::update(const reco::BeamSpot &beamspot,
                                         const float beam_weight,
                                         const bool fill_covariances) {
  double delta_z = 0;
  double delta_x = 0;
  double delta_y = 0;
  unsigned const int nt = trackinfo.size();
  unsigned const int nv = xv.size();
  if (fill_covariances) {
    V_vtx.clear();
  }

  // initial value for S, 0 or inverse of the beamspot covariance matrix
  Error3 S0;
  double c0 = 0, c1 = 0, c2 = 0;
  if (beam_weight > 0) {
    S0 = get_inverse_beam_covariance(beamspot);
    c0 = -(S0(0, 0) * beamspot.x0() + S0(0, 1) * beamspot.y0() + S0(0, 2) * beamspot.z0());
    c1 = -(S0(1, 0) * beamspot.x0() + S0(1, 1) * beamspot.y0() + S0(1, 2) * beamspot.z0());
    c2 = -(S0(2, 0) * beamspot.x0() + S0(2, 1) * beamspot.y0() + S0(2, 2) * beamspot.z0());
  }

  for (unsigned int k = 0; k < nv; k++) {
    Error3 S(S0);
    // sum track contributions
    double c[3] = {c0, c1, c2};
    for (unsigned int j = tkfirstv[k]; j < tkfirstv[k + 1]; j++) {
      const unsigned int i = tkmap[j];
      const auto w = tkweight[j];

      S += w * trackinfo[i].C;
      for (unsigned int l = 0; l < 3; l++) {
        c[l] += w * trackinfo[i].g[l];
      }
    }

#ifdef PVTX_DEBUG
    if ((fabs(S(1, 2) - S(2, 1)) > 1e-3) || (fabs(S(0, 2) - S(2, 0)) > 1e-3) || (fabs(S(0, 1) - S(1, 0)) > 1e-3) ||
        (S(0, 0) <= 0) || (S(0, 0) <= 0) || (S(0, 0) <= 0)) {
      edm::LogWarning("AMultiPrimaryVertexFitter") << "update()  bad S-matrix   S=" << std::endl << S << std::endl;
      edm::LogWarning("AMultiPrimaryVertexFitter") << "n-vertex = " << nv << "  n-track = " << nt << std::endl;
    }
#endif

    const auto xold = xv[k];
    const auto yold = yv[k];
    const auto zold = zv[k];

    if (S.Invert()) {
      xv[k] = -(S(0, 0) * c[0] + S(0, 1) * c[1] + S(0, 2) * c[2]);
      yv[k] = -(S(1, 0) * c[0] + S(1, 1) * c[1] + S(1, 2) * c[2]);
      zv[k] = -(S(2, 0) * c[0] + S(2, 1) * c[1] + S(2, 2) * c[2]);
      if (fill_covariances) {
        V_vtx.emplace_back(S);
      }
    } else {
      edm::LogWarning("AMultiPrimaryVertexFitter") << "update()   Matrix inversion failed" << S << std::endl;
    }

    if ((nt > 1) && (rho_vtx[k] > (1. / nt))) {
      delta_x = std::max(delta_x, std::abs(xv[k] - xold));
      delta_y = std::max(delta_y, std::abs(yv[k] - yold));
      delta_z = std::max(delta_z, std::abs(zv[k] - zold));
    }
  }  // vertex loop

  return std::max(delta_z, std::max(delta_x, delta_y));
}

AMultiPrimaryVertexFitter::Error3 AMultiPrimaryVertexFitter::get_inverse_beam_covariance(
    const reco::BeamSpot &beamspot) {
  auto SBeam = beamspot.rotatedCovariance3D();
  if (SBeam.Invert()) {
    return SBeam;
  } else {
    edm::LogWarning("AMultiPrimaryVertexFitter")
        << "Warning, beam-spot covariance matrix inversion failed " << std::endl;
    Error3 S0;
    S0(0, 0) = 1. / pow(beamspot.BeamWidthX(), 2);
    S0(1, 1) = 1. / pow(beamspot.BeamWidthY(), 2);
    S0(2, 2) = 1. / pow(beamspot.sigmaZ(), 2);
    return S0;
  }
}

std::vector<TransientVertex> AMultiPrimaryVertexFitter::fit(const std::vector<reco::TransientTrack> &tracks,
                                                            const std::vector<TransientVertex> &clusters,
                                                            const reco::BeamSpot &beamspot,
                                                            const bool useBeamConstraint) {
  const int max_iterations = 50;

  // initialize the vertices
  const unsigned int nv = clusters.size();
  xv.clear();
  xv.reserve(nv);
  yv.clear();
  yv.reserve(nv);
  zv.clear();
  zv.reserve(nv);
  rho_vtx.clear();
  rho_vtx.reserve(nv);
  tkfirstv.clear();
  tkfirstv.reserve(nv + 1);
  V_vtx.clear();
  V_vtx.reserve(nv);

  // seeds
  unsigned int sum_numtk = 0;
  for (auto &clu : clusters) {
    const double zclu = clu.position().z();
    xv.emplace_back(beamspot.x(zclu));
    yv.emplace_back(beamspot.y(zclu));
    zv.emplace_back(zclu);
    const unsigned int numtk = !clu.originalTracks().empty() ? clu.originalTracks().size() : 1;
    rho_vtx.emplace_back(numtk);
    sum_numtk += numtk;
  }
  // normalize the vertex multiplicites
  for (unsigned int k = 0; k < nv; k++) {
    rho_vtx[k] /= sum_numtk;
  }

  // use all tracks
  input_tracks = tracks;
  fill_trackinfo(tracks, beamspot);

  float beam_weight = useBeamConstraint ? 1. : 0.;

  double beta = 1.0;  // no annealing (don't undo what the clusterizer did)
  make_vtx_trk_map(beta);
  double delta = 0;
  unsigned int nit = 0;
  while ((nit == 0) || ((delta > 0.0001) && (nit < max_iterations))) {
    fill_weights(beta, beamspot);
    delta = update(beamspot, beam_weight, false);
    nit++;
  }
  if (nit >= max_iterations) {
    edm::LogWarning("AMultiPrimaryVertexFitter")
        << "iteration limit reached " << nit << "  last delta = " << delta << std::endl
        << " nv = " << nv << "    nt = " << tracks.size() << std::endl;
  }

  update(beamspot, beam_weight, true);  // fill the covariance matrices this time

  // track assignment : for each track identify the vertex that wants it most
  std::vector<unsigned int> kmaxweight(trackinfo.size(), nv);
  std::vector<double> maxweight(trackinfo.size(), -1.);
  for (unsigned int k = 0; k < nv; k++) {
    for (unsigned int j = tkfirstv[k]; j < tkfirstv[k + 1]; j++) {
      const unsigned int i = tkmap[j];
      if (tkweight[j] > maxweight[i]) {
        maxweight[i] = tkweight[j];
        kmaxweight[i] = k;
      }
    }
  }

  // fill the fit result into transient vertices
  std::vector<TransientVertex> pvs;

  for (unsigned int k = 0; k < xv.size(); k++) {
    const GlobalPoint pos(xv[k], yv[k], zv[k]);
    const GlobalError posError(
        V_vtx[k](0, 0), V_vtx[k](1, 0), V_vtx[k](1, 1), V_vtx[k](2, 0), V_vtx[k](2, 1), V_vtx[k](2, 2));
    float chi2 = 0.;
    float vtx_ndof = -3.;
    if (useBeamConstraint) {
      // add beam-spot chi**2 and degrees of freedom
      vtx_ndof = 0;
      const auto S = get_inverse_beam_covariance(beamspot);
      const double dx = xv[k] - beamspot.x0();
      const double dy = yv[k] - beamspot.y0();
      const double dz = zv[k] - beamspot.z0();
      chi2 = beam_weight * (S(0, 0) * dx * dx + S(1, 1) * dy * dy + 2 * S(0, 1) * dx * dy + S(2, 2) * dz * dz +
                            2 * S(0, 2) * dx * dz + 2 * S(1, 2) * dy * dz);
    }

    std::vector<reco::TransientTrack> vertex_tracks;
    TransientVertex::TransientTrackToFloatMap trkWeightMap;
    for (unsigned int j = tkfirstv[k]; j < tkfirstv[k + 1]; j++) {
      const unsigned int i = tkmap[j];
      const auto track_weight = tkweight[j];
      if ((track_weight >= min_trackweight_) && (kmaxweight[i] == k)) {
        vertex_tracks.emplace_back(input_tracks[i]);
        trkWeightMap[input_tracks[i]] = track_weight;
        vtx_ndof += 2 * track_weight;
        chi2 += track_weight * track_in_vertex_chsq(trackinfo[i], xv[k], yv[k], zv[k]);
      }
    }

    auto vtx = TransientVertex(pos, posError, vertex_tracks, chi2, vtx_ndof);
    vtx.weightMap(trkWeightMap);
    pvs.emplace_back(vtx);
  }

  return pvs;
}

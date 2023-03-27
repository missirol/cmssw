#include "RecoVertex/PrimaryVertexProducer/interface/MultiPrimaryVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

//#define DEBUG
#ifdef DEBUG
#define DEBUGLEVEL 0
#endif

/*
MultiPrimaryVertexFitter(const edm::ParameterSet &conf){
  chi2cutoff_ = conf.getUntrackedParameter<double>("chi2cutoff", 2.5);
}
*/

const bool verbose_ = false;

MultiPrimaryVertexFitter::MultiPrimaryVertexFitter(const double chi2cutoff, const double mintrkweight):chi2_cutoff_(chi2cutoff), min_trackweight_(mintrkweight){
  if(verbose_){
    std::cout << "instantiating MPV with   chi2_cutoff =" <<  chi2_cutoff_ << "  and min_trackweight=" << min_trackweight_ << std::endl;
  }
}


void MultiPrimaryVertexFitter::fill_trackinfo(const std::vector<reco::TransientTrack>  & tracks, const reco::BeamSpot &beamSpot){
  /* fill track information used during fits into arrays, parallell to the list of input tracks */
  double zB = beamSpot.z0();
  double xB = beamSpot.x0();
  double yB = beamSpot.y0();

  trackinfo.clear();
  trackinfo.reserve(tracks.size());

  for(auto & trk : tracks){
    TrackInfo ti;
    const auto & p0 = trk.track().parameters();
    const auto & cov = trk.track().covariance();
    const auto phi = p0[reco::TrackBase::i_phi];
    const auto lambda = p0[reco::TrackBase::i_lambda];
    const auto dsz = p0[reco::TrackBase::i_dsz];
    const auto dxy = p0[reco::TrackBase::i_dxy];
    auto const cos_phi = cos(phi);
    auto const sin_phi = sin(phi);
    auto cos_lambda = cos(lambda);
    auto sin_lambda = sin(lambda);

    // the distances of closest approach of a track to an arbitrary point (x,y,z)
    // in the straight track approximation are
    // F1() = dxy() = dxy + x * sin_phi - y cos_phi
    // F2() = dsz() = dsz  - z * cos_lambda + (x cos_phi + y sin_phi) * sin_lambda
    // it is assumed, that x,y are reasonably close to the beam-spot (x=xb+x', y=yb+y')
    // F (x) = F(xB + x') = b + a x'
    // error propagation:
    // J12 = dF1 / dpar(2) ...,  par(1)=lambda, par(2)=phi, par(3)=dxy, par(4)=dsz
    double J12 = xB * cos_phi + yB * sin_phi;
    double J13 = 1.;
    double J21 = zB * sin_lambda + (xB * cos_phi + yB * sin_phi) * cos_lambda;
    double J22 = (xB * sin_phi + yB * cos_phi) * sin_lambda;
    double J24 = 1.;
    // U_{a,b} = cov<F_a, F_b> = sum_{u,v} J_{au} J_{bv} V_{uv}, // J11=J23=J14=0
    double U11 = J12 * J12 * cov(2,2) + J13 * J13 * cov(3,3) + 2 * J12 * J13 * cov(2,3); 
    double U22 = J21 * J21 * cov(1,1) + J22 * J22 * cov(2,2) + J24 * J24 * cov(4,4)
      + 2 * (J21 * J22 * cov(1,2) + J21 * J24 * cov(1,4) + J22 * J24 * cov(2,4));
    double U12 =  J12 * J22 * cov(2,2)                                                    
      + J12 * J21 * cov(1,2) + J12 * J24 * cov(2,4)
      + J13 * J21 * cov(1,3) + J13 * J22 * cov(3,2) + J13 * J24 * cov(3,4);
    // S = U^{-1}
    double DetU = U11 * U22 - U12 * U12;
    if(fabs(DetU)<1.e-16){
      std::cout << "Warning, det(U) almost vanishes : " << DetU << " !! This should not happen!" << std::endl;
      ti.S11 = 0;
      ti.S22 = 0;
      ti.S12 = 0;
    }else{
      ti.S11 =   U22 / DetU;
      ti.S22 =   U11 / DetU;
      ti.S12 = - U12 / DetU;
    }
#ifdef DEBUG    
    assert((ti.S11 >=0)&&"S11 positivity test");
    assert((ti.S22 >=0)&&"S22 positivity test");
#endif    

    ti.b1 = dxy + sin_phi * xB - cos_phi * yB;
    ti.a1[0] =   sin_phi;
    ti.a1[1] = - cos_phi;
    ti.a1[2] = 0;
    ti.b2 = dsz - zB * cos_lambda + (cos_phi * xB + sin_phi * yB) * sin_lambda;
    ti.a2[0] = cos_phi * sin_lambda;
    ti.a2[1] = sin_phi * sin_lambda;
    ti.a2[2] = -cos_lambda;

    // c = b^T S A  [2^T x (2,2) x (2,3) -> 3]     = 2^T x (2,3)
    // C = A^T S A  [(3,2) x (2,2) x (2,3) -> 3x3] = (3,2) x ( 2,3) (symmetric)
    for(int k = 0; k < 3; k++){
      double SA1k = (ti.S11 * ti.a1[k] + ti.S12 * ti.a2[k]);
      double SA2k = (ti.S12 * ti.a1[k] + ti.S22 * ti.a2[k]);
      ti.c[k] = ti.b1 * SA1k + ti.b2 * SA2k;
      for(int l = 0; l < 3; l++){
        ti.C(l,k) = ti.a1[l] * SA1k + ti.a2[l] * SA2k;
      }
    }
    ti.d = ti.b1 * ti.b1 * ti.S11 + ti.b2 * ti.b2 * ti.S22 + 2 * ti.b1 * ti.b2 * ti.S12;

    assert((fabs(ti.C(1,2)-ti.C(2,1))<1e-3) && " C12 symmetry test");
    assert((fabs(ti.C(0,2)-ti.C(2,0))<1e-3) && " C02 symmetry test");
    assert((fabs(ti.C(0,1)-ti.C(1,0))<1e-3) && " C01 symmetry test");
    assert((ti.C(0,0)>0) && " C00  positivity test");
    assert((ti.C(1,1)>0) && " C11  positivity test");
    assert((ti.C(2,2)>0) && " C22  positivity test");
    
    ti.z = trk.track().vz();
    ti.x = trk.track().vx(); // for debugging
    ti.y = trk.track().vy(); // for debugging
    ti.ipsig = trk.stateAtBeamLine().transverseImpactParameter().significance();// for debugging
    ti.odz2 = 1. / pow(trk.track().dzError(), 2);

    trackinfo.push_back(ti);
  }
}



void MultiPrimaryVertexFitter::fill_weights(double beta,  const reco::BeamSpot & beamspot, const double Zcutoff){
  unsigned const int nt = trackinfo.size();
  unsigned const int nv = xv.size();

  const auto xb = beamspot.x0();
  const auto yb = beamspot.y0();
  const auto zb = beamspot.z0();
  double rho_sum = 0;
  double Z_noise = 0;
  double rho_noise = 0;
  for(unsigned int k = 0; k < nv; k++){
    chi2_vtx[k] = 0;
    d2Fperp_vtx[k] = 0;
    rho_sum += rho_vtx[k];
  }

  if(rho_sum < 0.99999){
    rho_noise = 1. - rho_sum;
    Z_noise = rho_noise * exp(-0.5 * beta * chi2_cutoff_ * 2 - log(nv));
  }else{
    rho_noise = 0;
    Z_noise = 0;
  }
  
  std::vector<double> Z(nt, 0.);
  std::vector<double> ccache(nv, 0.);
  std::vector<double> ecache(nv, 0.);
  std::vector<double> rhocache(nv, 0.);
  double rhoXcache = 0; // debugging
  std::vector<double> d2Fperpcache(nv, 0.);

  // temp debug
  std::vector<double> temp_chsqmin(nt, 1e10);
  std::vector<double> temp_chsqzmin(nt, 1e10);
  

  for(unsigned int i=0; i < nt; i++){
    auto & ti = trackinfo[i];
    ti.weight.clear();
    ti.weight.reserve(nv);

    // evaluate and cache track-vertex assignment chi**2 for all clusters and sum up Z[i]
    for(unsigned int k = 0; k < nv; k++){
      ti.weight.emplace_back(0);
      double F1 = ti.b1 + (xv[k] - xb) * ti.a1[0] + (yv[k] - yb) * ti.a1[1] ; // a1[2]==0
      double F2 = ti.b2 + (xv[k] - xb) * ti.a2[0] + (yv[k] - yb) * ti.a2[1] +  (zv[k] - zb) * ti.a2[2];
      double chsq = F1 * F1 * ti.S11 + F2 * F2 * ti.S22 + 2. * F1 * F2 * ti.S12;
      assert((chsq >= 0) && " negative chi**2");
      ccache[k] = chsq;
      // notes : 
      //                           F2 * F2 * ti.S22                       approx=      pow(ti.z-zv[k], 2) * ti.odz2 
      // hence           (1. - beta * F2 * F2 * ti.S22)  * ti.odz2        approx=     (1. - beta * pow(ti.z-zv[k], 2)* ti.odz2) * ti.odz2
      // but   different from  (1. - beta * F2 * F2 * ti.S22)  * ti.S22
      d2Fperpcache[k] = (1. - beta * pow(ti.z-zv[k], 2)* ti.odz2) * ti.odz2; // lacks a factor p_{ik} -> next loop 

      if(chsq < temp_chsqmin[i]) temp_chsqmin[i] = chsq;
      double chsqz = F2 * F2 * ti.S22;
      if(chsqz < temp_chsqzmin[i]) temp_chsqzmin[i] = chsqz;
      
      if((beta * chsq) < chi2_cutoff_){
	double arg = 0.5 * beta * chsq;
	ecache[k] = exp(-arg);
        Z[i]  += rho_vtx[k] * ecache[k];
      }else{
	ecache[k] = 0.;
	d2Fperpcache[k] = 0;
      }
      
    }// end of the 1st cluster loop for track i


    // now we have the partition function, Z_i and can evaluate assignment probabilities, loop again over clusters (for track i)
    if (Z[i] > 0){
      rhoXcache += Z_noise / (Z[i] + Z_noise);
      for(unsigned int k = 0; k < nv; k++){
	double weight = rho_vtx[k] * ecache[k] / (Z[i] + Z_noise);
	ti.weight[k] = weight;
	chi2_vtx[k] += ccache[k] * weight;
	rhocache[k] += weight;
	d2Fperp_vtx[k] += weight * d2Fperpcache[k];
      }
    }else{
      for(unsigned int k = 0; k < nv; k++){
	ti.weight[k] = 0.;// redundant
      }
      rhoXcache += 1.;
      if(verbose_){
	std::cout <<" updateWeights : Warning ! track " << i 
		  << " ztrk=" <<   std::setw(10) << std::fixed << std::setprecision(4) << ti.z  
		  << " Z=0: " << Z[i] 
                  <<   std::setw(5) << std::fixed  <<  std::setprecision(1) << "  chi**2 min = " << temp_chsqmin[i] 
		  <<  "   chi2zmin = " << temp_chsqzmin[i] << "   chi2_cutoff_=" << chi2_cutoff_ << "  ipsig=" << ti.ipsig << std::endl;
      }
    }

  }//track loop



  double sum_rho = 0.;
  for(unsigned int k=0; k < nv; k++){
    rho_vtx[k] = rhocache[k] / nt;
    sum_rho += rho_vtx[k];
  }
  // just a test, is the noise cluster behaviour consistent?
  if( std::abs(1.-sum_rho - rhoXcache / nt) > 0.01){
    std::cout << " fill_weights sum_rho = " << sum_rho  << "  1-sum_rho = " <<  1.-sum_rho <<   " rhoX=" << rhoXcache / nt  << std::endl;
  }

}


/**********************************************************************************/
/*
void MultiPrimaryVertexFitter::test_chisquared(){
      // test chi**2 = x'T C x' + 2 cTx' + d  == FT S F
      const double tx = xv[k] - xb;
      const double ty = yv[k] - yb;
      const double tz = zv[k] - zb;
      // makes chi**2 vanish :
      //const double tx = ti.x -xb;
      //const double ty = ti.y -yb;
      //const double tz = ti.z -zb;
      double dchsq1 =
            pow(tx, 2) * ti.C(0,0)
	+   pow(ty, 2) * ti.C(1,1)
	+   pow(tz, 2) * ti.C(2,2)
	+   tx * ty * (ti.C(0,1) + ti.C(1,0))
	+   tx * tz * (ti.C(0,2) + ti.C(2,0))
	+   ty * tz * (ti.C(1,2) + ti.C(2,1));
      double dchsq2 = 
	2* (ti.c[0] * tx + ti.c[1] * ty + ti.c[2] * tz);
      
      double tF1 = ti.b1 + tx * ti.a1[0] + ty * ti.a1[1] ; 
      double tF2 = ti.b2 + tx * ti.a2[0] + ty * ti.a2[1] +  tz * ti.a2[2];
      double tchsq = tF1 * tF1 * ti.S11 + tF2 * tF2 * ti.S22 + 2. * tF1 * tF2 * ti.S12 ;

      std::cout << "YYYYY  xychsq = " <<   tchsq << "   dchsq = "<< dchsq1+dchsq2+ti.d << "   dchsq1=" << dchsq1  << "   dchsq2=" << dchsq2  << "   dchsq3=" << ti.d  <<  "  zb=" << zb << std::endl;
      std::cout << "const  " << ti.d << " " << ti.b1*ti.b1*ti.S11 + ti.b2*ti.b2*ti.S22  + 2*ti.b1*ti.b2*ti.S12 << std::endl;
      std::cout << "tx     " << ti.c[0] << " " << ti.b1*ti.a1[0]*ti.S11 + ti.b2*ti.a2[0]*ti.S22 + ti.b1*ti.a2[0]*ti.S12 + ti.a1[0]*ti.b2*ti.S12 << std::endl;
      std::cout << "ty     " << ti.c[1] << " " << ti.b1*ti.a1[1]*ti.S11 + ti.b2*ti.a2[1]*ti.S22 + ti.b1*ti.a2[1]*ti.S12 + ti.a1[1]*ti.b2*ti.S12 << std::endl;
      std::cout << "tz     " << ti.c[2] << " " << ti.b2*ti.a2[2]*ti.S22 + ti.b1*ti.a2[2]*ti.S12 << std::endl;
      std::cout << "txtx   " << ti.C(0,0) << " " << ti.a1[0]*ti.a1[0]*ti.S11 + ti.a2[0]*ti.a2[0]*ti.S22 + 2*ti.a1[0]*ti.a2[0]*ti.S12 << std::endl;
      std::cout << "tyty   " << ti.C(1,1) << " " << ti.a1[1]*ti.a1[1]*ti.S11 + ti.a2[1]*ti.a2[1]*ti.S22 + 2*ti.a1[1]*ti.a2[1]*ti.S12 << std::endl;
      std::cout << "tztz   " << ti.C(2,2) << " " << ti.a2[2]*ti.a2[2]*ti.S22 << std::endl;
      std::cout << "txty   " << ti.C(0,1) << "," << ti.C(1,0) << " " << ti.a1[0]*ti.a1[1]*ti.S11 + ti.a2[0]*ti.a2[1]*ti.S22 + ti.a1[0]*ti.a2[1]*ti.S12 + ti.a1[1]*ti.a2[0]*ti.S12 << std::endl;
      std::cout << "txtz   " << ti.C(0,2) << "," << ti.C(2,0) << " " << ti.a2[0]*ti.a2[2]*ti.S22 + ti.a1[0]*ti.a2[2]*ti.S12 << std::endl;
      std::cout << "tytz   " << ti.C(1,2) << "," << ti.C(1,2 )<< " " << ti.a2[1]*ti.a2[2]*ti.S22 + ti.a1[1]*ti.a2[2]*ti.S12 << std::endl;
 }
*/
/**********************************************************************************/



double MultiPrimaryVertexFitter::single_fit(const reco::BeamSpot & beamspot, const float beam_weight, const bool fill_covariances){
  double delta_z = 0;
  double delta_x = 0;
  double delta_y = 0;
  unsigned const int nt = trackinfo.size();
  unsigned const int nv = xv.size();
#ifdef DEBUG
  assert((nv== yv.size()) && "nv y mismatch");
  assert((nv== zv.size()) && "nv z mismatch");
  assert((nv== rho_vtx.size()) && "nv rho mismatch");
  assert((nv== chi2_vtx.size()) && "nv chi2 mismatch");
#endif
  if(fill_covariances) {
    V_vtx.clear();
  }

  // initial value for S, 0 or inverse of the beamspot covariance matrix
  Error3 S0;
  if(beam_weight > 0){
    auto SBeam = beamspot.rotatedCovariance3D();
    if (SBeam.Invert()){
      S0 = beam_weight * SBeam;
    }else{
      std::cout << "Warning, beam-spot covariance matrix inversion failed " << std::endl;
      S0(0,0)= beam_weight / pow(beamspot.BeamWidthX(), 2);
      S0(1,1)= beam_weight / pow(beamspot.BeamWidthY(), 2);
      S0(2,2)= beam_weight / pow(beamspot.sigmaZ(), 2);
    }
  }


  for(unsigned int k = 0; k < nv; k++){

    Error3 S(S0);
    // sum track contributions
    double c[3] = {0};
    for(unsigned int i = 0; i < nt; i++){
      const double w = trackinfo[i].weight[k];
      S += w * trackinfo[i].C;
      for(unsigned int j = 0; j < 3; j++){
	c[j] += w * trackinfo[i].c[j];
      }
    }

#ifdef DEBUG    
    if((fabs(S(1,2)-S(2,1))>1e-3) || (fabs(S(0,2)-S(2,0))>1e-3) || (fabs(S(0,1)-S(1,0))>1e-3) || (S(0,0)<=0)  || (S(0,0)<=0)  || (S(0,0)<=0)){
      std::cout << "MultiPrimaryVertexFitter::single_fit  bad S-matrix   S=" << std::endl << S << std::endl;
      std::cout << "n-vertex = " << nv << "  n-track = " << nt << std::endl;
    }
#endif
    
    const auto xold = xv[k];
    const auto yold = yv[k];
    const auto zold = zv[k];

    if(S.Invert()){
      xv[k] = beamspot.x0() - ( S(0,0) * c[0] + S(0,1) * c[1] + S(0,2) * c[2] );
      yv[k] = beamspot.y0() - ( S(1,0) * c[0] + S(1,1) * c[1] + S(1,2) * c[2] );
      zv[k] = beamspot.z0() - ( S(2,0) * c[0] + S(2,1) * c[1] + S(2,2) * c[2] );
      if(fill_covariances){
	V_vtx.emplace_back(S);
      }
    }else{
      std:: cout <<"MultiPrimaryVertexFitter::single_fit   Matrix inversion failed" << S << std::endl;
    }

    if((rho_vtx[k] * nt) > 1.0){
      delta_x = std::max(delta_x, std::abs(xv[k] - xold));
      delta_y = std::max(delta_y, std::abs(yv[k] - yold));
      delta_z = std::max(delta_z, std::abs(zv[k] - zold));
    }
    
  }// vertex loop

  return std::max(delta_z, std::max(delta_x,delta_y));
}


void MultiPrimaryVertexFitter::dump(const std::string & comment, const::std::vector<TransientVertex> & clusters, 
				    const std::vector<reco::TransientTrack> & tracks,
				    const reco::BeamSpot & beamspot,
				    const double zmin, const double zmax, const unsigned int nit=666){
  double sumrho = 0;
  std::cout << "dump ----------------- " << comment << "   nv=" << zv.size() << "   nt=" << trackinfo.size();
  if(nit == 666) { 
    std::cout << std::endl;
  }else{
    std::cout << "  nit=" << nit << std::endl;
  }
  const auto xb = beamspot.x0();
  const auto yb = beamspot.y0();
  const auto zb = beamspot.z0();
  
  for(unsigned int k=0; k < zv.size(); k++){
    sumrho+=rho_vtx[k];
    if ( (zv[k] > zmin) && (zv[k] < zmax) ){
      std::cout <<  std::setw(10) << std::setprecision(4) << std::fixed << clusters[k].position().z() 
		<< "  ->  " 
		<< std::setw(10) << std::setprecision(4) << std::fixed << zv[k]  
		<< "  x,y= " 
		<< std::setw(10) << std::setprecision(4) << std::fixed << xv[k] << "," 
		<< std::setw(10) << std::setprecision(4) << std::fixed << yv[k] 
		<< "     rho*nt= "
		<< std::setw(10) << std::setprecision(3) << std::fixed << rho_vtx[k]  * trackinfo.size()
		<< "     rho*nv= "
		<< std::setw(10) << std::setprecision(3) << std::fixed << rho_vtx[k]  * zv.size()
		<< "     nt(cluster) " 
		<< std::setw(10) << std::fixed << clusters[k].originalTracks().size();
      std::cout << "     d2Fperp = " << std::setw(12) << std::scientific << d2Fperp_vtx[k];
      if(d2Fperp_vtx[k] < 0){
	std::cout << " MERGER?";
      }else{
	std::cout << "        ";
      }

      std::cout << std::endl;

      // dump track details on selected vertices
      //if((k==555555) ||((clusters[k].originalTracks().size() > 5) && ((fabs(xv[k]) > 0.025) || (fabs(yv[k]) > 0.025)))){
      //	std::cout << " far-out vertex" << k << std::endl;
      if( k < 5){
	for(unsigned int i=0; i<trackinfo.size(); i++){
	  // is it an original part of the cluster?
	  bool in_original_tracks = false;
	  for(const auto & tk : clusters[k].originalTracks()){
	      if (((std::abs(trackinfo[i].z - tk.track().vz())<0.00001)) && ((std::abs(trackinfo[i].x - tk.track().vx())<0.00001)) ){
		in_original_tracks = true;
	      }
	  }
	  
	  auto & ti = trackinfo[i];
	  if((ti.weight[k] < 0.01) && (!in_original_tracks) )continue;
	  if(in_original_tracks){
	    std::cout << "  * " ;
	  }else{
	    std::cout << "    " ;
	  }

	  double F1 = ti.b1 + (xv[k] - xb) * ti.a1[0] + (yv[k] - yb) * ti.a1[1] ; // a1[2]==0
	  double F2 = ti.b2 + (xv[k] - xb) * ti.a2[0] + (yv[k] - yb) * ti.a2[1] +  (zv[k] - zb) * ti.a2[2];
	  double chsqv = F1 * F1 * ti.S11 + F2 * F2 * ti.S22 + 2. * F1 * F2 * ti.S12;
	  F1 = ti.b1;
	  F2 = ti.b2 +  (zv[k] - zb) * ti.a2[2];
	  double chsq0 = F1 * F1 * ti.S11 + F2 * F2 * ti.S22 + 2. * F1 * F2 * ti.S12;

	  std::cout << std::setw(5) << std::fixed  <<  std::setprecision(3) <<ti.weight[k] 
		    << std::setw(8) << std::fixed  << std::setprecision(4) << ti.x
		    << std::setw(8) << std::fixed  << std::setprecision(4) << ti.y
		    << std::setw(10) << std::fixed << std::setprecision(4) << ti.z;
	    std::cout << "  ez= " << std::setw(7) << std::setprecision(4) << tracks[i].track().dzError()
		      << "  exy= "  << std::setw(7) << std::setprecision(4) << tracks[i].track().dxyError()
		      << "  phi= "  << std::setw(7) << std::setprecision(4) << tracks[i].track().phi()
		      << "  eta= "  << std::setw(7) << std::setprecision(4) << tracks[i].track().eta()
		      << "  pt= "  << std::setw(7) << std::setprecision(4) << tracks[i].track().pt()
		      << "  sxy(b)= " << std::setw(7) << std::setprecision(4) << tracks[i].track().dxy(beamspot)/ tracks[i].track().dxyError()
		      << "  sxy(v)= " << std::setw(7) << std::setprecision(4) << tracks[i].track().dxy(math::XYZPoint(xv[k],yv[k],zv[k]))/ tracks[i].track().dxyError()
		      << "  chsqv= " <<  std::setw(5) << std::setprecision(1) << chsqv
		      << "  chsq0= " <<  std::setw(5) << std::setprecision(1) << chsq0
		      << "  wdchsq= " <<  std::setw(5) << std::setprecision(1) << (chsq0-chsqv) * ti.weight[k]
	      ;
	  std::cout << std::endl;
	}

      }


      // and now dump a chi**2 map
      if( (hdump < 0) && (
	 ((clusters[k].originalTracks().size() > 5) && ((fabs(xv[k]) > 0.025) || (fabs(yv[k]) > 0.025)) && (fabs(zv[k]-0.11)<0.2))
	 || (k==55))
	 ){
	hdump++;
	double z = zv[k];
	std::cout << " dumping  the chi**2 map for z= " << z << "   dump " << hdump << std::endl;
	for(int row=500; row > -501; row -= 10){
	  double y = row * 1e-4;
	  for(int col=-500; col<501; col += 10){
	    double x = col * 1e-4;
	    double sumchsq = 0;
	    
	    for(unsigned int i=0; i < trackinfo.size(); i++){
	      const auto ti = trackinfo[i];
	      if(ti.weight[k] < 0.01) continue;
	      double F1 = ti.b1 + (x - xb) * ti.a1[0] + (y - yb) * ti.a1[1] ; // a1[2]==0
	      double F2 = ti.b2 + (x - xb) * ti.a2[0] + (y - yb) * ti.a2[1] +  (z - zb) * ti.a2[2];
	      double chsq = F1 * F1 * ti.S11 + F2 * F2 * ti.S22 + 2. * F1 * F2 * ti.S12;
	      sumchsq += ti.weight[k] * chsq;
	    }
	    
	    std::cout << "HH" << hdump << " " 
		      << std::fixed << std::setw(5) << row << std::setw(5) << col 
		      << std::setw(10) << std::setprecision(4) << x 
		      << std::setw(10) << std::setprecision(4) << y
		      << std::scientific << std::setw(12) << sumchsq 
		      << std::endl;
	  }
	}
      }
    }
  }
  if(zv.size() > 1){
    std::cout << " sum rho = " << sumrho << std::endl;
  }

}




void MultiPrimaryVertexFitter::clean(const std::vector<reco::TransientTrack> & tracks, const std::vector<TransientVertex> & clusters){
// clean up the cluster list (dummy for now)
  if(clusters.size()<2) return;
  for(unsigned int k1 = 0; k1 < clusters.size()-1; k1++){
    for(unsigned int k2 = k1+1; k2 < clusters.size(); k2++){
      if(std::abs(zv[k1]-zv[k2]) > 0.0100) continue;
      // test overlap
      double m1=0, m2=0, m11=0, m22=0;
      for(unsigned int i=0; i<trackinfo.size(); i++){
	double w1 = trackinfo[i].weight[k1];
	double w2 = trackinfo[i].weight[k2];
	m1 += w1;
	m2 += w2;
	if(w1>w2){
	  m11 += w1;
	}
	if(w2>w1){
	  m22 += w2;
	}
      }
    }
  }
}



std::vector<TransientVertex> MultiPrimaryVertexFitter::fit(const std::vector<reco::TransientTrack> & tracks, const std::vector<TransientVertex> & clusters, const reco::BeamSpot & beamspot, const bool useBeamConstraint){

  // normalized number of tracks for initial guess of rho
  unsigned int numtk = 0;
  for(auto & clu : clusters){
    numtk += clu.originalTracks().size();
  }
  hdump =0;

  // initialize the vertices
  xv.clear();  xv.reserve(clusters.size());
  yv.clear();  yv.reserve(clusters.size());
  zv.clear();  zv.reserve(clusters.size());
  rho_vtx.clear(); rho_vtx.reserve(clusters.size());
  V_vtx.clear(); V_vtx.reserve(clusters.size());
  chi2_vtx.clear(); chi2_vtx.reserve(clusters.size());
  d2Fperp_vtx.clear(); d2Fperp_vtx.reserve(clusters.size());

  // seeds
  for(auto & clu : clusters){
    const double zclu = clu.position().z();
    xv.emplace_back(beamspot.x(zclu));
    yv.emplace_back(beamspot.y(zclu));
    zv.emplace_back(zclu);
    chi2_vtx.emplace_back(0.);
    d2Fperp_vtx.emplace_back(0);
    rho_vtx.emplace_back(double(clu.originalTracks().size()) / numtk);
  }

  // use all tracks
  input_tracks = tracks;
  fill_trackinfo(tracks, beamspot);


  // no annealing for now
  float beta = 1.0;
  double delta = 0;
  double Zcutoff = 0.;
  float beam_weight = 1.; 

  unsigned int nit = 0;
  while((nit==0) || ((delta > 0.0001) && (nit < 50))){
    fill_weights(beta, beamspot);
    delta = single_fit(beamspot, beam_weight, true);
    nit++;
  }

  if(verbose_){
    if(clusters.size() > 1){
      std::cout << "MultiPrimaryVertexFitter::fit  converged in " << nit << " iterations   delta=" << delta
		<< "  clusters.size=" << clusters.size() << std::endl;
    }
    dump("final", clusters, tracks, beamspot, -20., 20., nit);
  }


  // track assignment : for each track: identify the vertex that wants it most
  std::vector<unsigned int> kmaxweight(trackinfo.size(), xv.size());
  for(unsigned int i = 0; i< trackinfo.size(); i++){
    double maxweight = -1.;
    for(unsigned int k = 0; k < zv.size(); k++){
      if(trackinfo[i].weight[k] > maxweight){
	maxweight = trackinfo[i].weight[k];
	kmaxweight[i] = k;
      }
    }
  }
#ifdef DEBUG
  for(unsigned int i = 0; i< trackinfo.size(); i++){
    if(maxweight > 1.1){
      assert(maxweight == trackinfo[i].weight[kmaxweight[i]]);
      const auto xb = beamspot.x0();
      const auto yb = beamspot.y0();
      const auto zb = beamspot.z0();
      const auto k = kmaxweight[i];
      auto & ti = trackinfo[i];
      double F1 = ti.b1 + (xv[k] - xb) * ti.a1[0] + (yv[k] - yb) * ti.a1[1] ; // a1[2]==0
      double F2 = ti.b2 + (xv[k] - xb) * ti.a2[0] + (yv[k] - yb) * ti.a2[1] +  (zv[k] - zb) * ti.a2[2];
      double chsq = F1 * F1 * ti.S11 + F2 * F2 * ti.S22 + 2. * F1 * F2 * ti.S12;

      std::cout << "track " << i <<  " vertex " << kmaxweight[i] << "  weight " << trackinfo[i].weight[kmaxweight[i]] 
		<< " rho = " << rho_vtx[kmaxweight[i]] 
		<< " maxweight = " << rho_vtx[kmaxweight[i]] / (Zcutoff/zv.size() +  rho_vtx[kmaxweight[i]])
		<< " Zcutoff_eff = " << Zcutoff/(zv.size()*rho_vtx[kmaxweight[i]])
		<< " Zcutoff = " << Zcutoff
		<< " F1 = " << F1 << " (" << F1*sqrt(ti.S11) <<")"
		<< " F2 = " << F2 << " (" << F2*sqrt(ti.S22) <<")" 
		<< " chsq = " << chsq
		<< " corr" << ti.S12/sqrt(ti.S11*ti.S22)
		<< std::endl;
    }
  }
#endif


  // fill the fit result into transient vetices
  std::vector<TransientVertex> pvs;

  for(unsigned int k = 0; k < xv.size(); k++){
    const GlobalPoint pos(xv[k], yv[k], zv[k]);
    const GlobalError posError(V_vtx[k](0,0), V_vtx[k](1, 0), V_vtx[k](1, 1), V_vtx[k](2, 0), V_vtx[k](2, 1), V_vtx[k](2,2));
    const float chi2 = chi2_vtx[k];
    float vtx_ndof = -3.;
    if(useBeamConstraint){
      vtx_ndof = 0.;
    }
    
    std::vector<reco::TransientTrack> vertex_tracks;
    TransientVertex::TransientTrackToFloatMap trkWeightMap;
    for(unsigned int i = 0; i < trackinfo.size(); i++){
      const auto track_weight = trackinfo[i].weight[k];
      if((track_weight > min_trackweight_) && (kmaxweight[i] == k)) {
	vertex_tracks.emplace_back(input_tracks[i]);
	trkWeightMap[input_tracks[i]] = track_weight;
	vtx_ndof += 2 * track_weight;
      }
    }
    
    auto vtx = TransientVertex(pos, posError, vertex_tracks, chi2, vtx_ndof);
    vtx.weightMap(trkWeightMap);
    pvs.emplace_back(vtx);
  }


  return pvs;
}


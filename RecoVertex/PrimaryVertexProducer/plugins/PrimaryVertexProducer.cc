// -*- C++ -*-
//
// Package:    PrimaryVertexProducer
// Class:      PrimaryVertexProducer
//
/**\class PrimaryVertexProducer PrimaryVertexProducer.cc RecoVertex/PrimaryVertexProducer/src/PrimaryVertexProducer.cc

 Description: steers tracker primary vertex reconstruction and storage

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Pascal Vanlaer
//         Created:  Tue Feb 28 11:06:34 CET 2006
//
//
#include <memory>
#include <algorithm>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZT_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/HITrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/SequentialPrimaryVertexFitterAdapter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/AdaptiveChisquarePrimaryVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/MultiPrimaryVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmFromTracksPID.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexTimeAlgorithmLegacy4D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/WeightedMeanFitter.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

class PrimaryVertexProducer : public edm::stream::EDProducer<> {
public:
  PrimaryVertexProducer(const edm::ParameterSet&);
  ~PrimaryVertexProducer() override = default;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  edm::EDGetTokenT<reco::BeamSpot> const bsToken;
  edm::EDGetTokenT<reco::TrackCollection> const trkToken;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> const theTTBToken;
  bool const fVerbose;

  std::unique_ptr<TrackFilterForPVFindingBase> theTrackFilter;
  std::unique_ptr<TrackClusterizerInZ> theTrackClusterizer;

  // vtx fitting algorithms
  struct algo {
    std::string label;
    bool useBeamConstraint;
    double minNdof;
    bool is4D;
    std::unique_ptr<VertexCompatibleWithBeam> pv_selector;
    std::unique_ptr<PrimaryVertexFitterBase> pv_fitter;
    std::unique_ptr<VertexTimeAlgorithmBase> pv_time_estimator;
  };

  std::vector<algo> algorithms;

  bool fRecoveryIteration;
  edm::EDGetTokenT<reco::VertexCollection> recoveryVtxToken;

  edm::EDGetTokenT<edm::ValueMap<float> > trkTimesToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeResosToken;

  bool useTransientTrackTime;
};

PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf) :
bsToken{consumes(conf.getParameter<edm::InputTag>("beamSpotLabel"))},
trkToken{consumes(conf.getParameter<edm::InputTag>("TrackLabel"))},
theTTBToken{esConsumes(conf.getParameter<edm::ESInputTag>("transientTrackBuilder"))},
fVerbose{conf.getUntrackedParameter<bool>("verbose")} {

  useTransientTrackTime = false;

  // select and configure the track selection
  std::string trackSelectionAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = std::make_unique<TrackFilterForPVFinding>(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = std::make_unique<HITrackFilterForPVFinding>(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }

  // select and configure the track clusterizer
  std::string clusteringAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm == "gap") {
    theTrackClusterizer = std::make_unique<GapClusterizerInZ>(conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  } else if (clusteringAlgorithm == "DA") {
    theTrackClusterizer = std::make_unique<DAClusterizerInZ>(conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }
  // provide the vectorized version of the clusterizer, if supported by the build
  else if (clusteringAlgorithm == "DA_vect") {
    theTrackClusterizer = std::make_unique<DAClusterizerInZ_vect>(conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }
  // 2D (z,t) Deterministic Annealing, vectorized
  else if (clusteringAlgorithm == "DA2D_vect") {
    theTrackClusterizer = std::make_unique<DAClusterizerInZT_vect>(conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
    useTransientTrackTime = true;
  }
  // unknown clustering algorithm, throw exception
  else {
    throw VertexException("PrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
  }

  // select and configure the vertex fitters
  auto const& vertexCollections = conf.getParameter<std::vector<edm::ParameterSet> >("vertexCollections");

  algorithms.clear();
  algorithms.reserve(vertexCollections.size());

  for (auto const& algoconf : vertexCollections) {
    algorithms.emplace_back();
    auto& algorithm = algorithms.back();

    // configure the fitter and selector
    auto const fitterAlgorithm = algoconf.getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.pv_fitter = std::make_unique<SequentialPrimaryVertexFitterAdapter>(std::make_unique<KalmanVertexFitter>());
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      algorithm.pv_fitter = std::make_unique<SequentialPrimaryVertexFitterAdapter>(std::make_unique<AdaptiveVertexFitter>(GeometricAnnealing(algoconf.getParameter<double>("chi2cutoff"))));
    } else if (fitterAlgorithm == "AdaptiveChisquareVertexFitter") {
      algorithm.pv_fitter = std::make_unique<AdaptiveChisquarePrimaryVertexFitter>(algoconf.getParameter<double>("chi2cutoff"), 0.);
    } else if (fitterAlgorithm == "MultiPrimaryVertexFitter") {
      algorithm.pv_fitter = std::make_unique<MultiPrimaryVertexFitter>(algoconf.getParameter<double>("chi2cutoff"), algoconf.getParameter<double>("mintrkweight"));
    } else if (fitterAlgorithm == "WeightedMeanFitter") {
      algorithm.pv_fitter = std::make_unique<WeightedMeanPrimaryVertexEstimator>();
    } else if (not fitterAlgorithm.empty()) {
      throw VertexException("PrimaryVertexProducer: unknown algorithm: " + fitterAlgorithm);
    }
    algorithm.label = algoconf.getParameter<std::string>("label");
    algorithm.minNdof = algoconf.getParameter<double>("minNdof");
    algorithm.useBeamConstraint = algoconf.getParameter<bool>("useBeamConstraint");
    algorithm.pv_selector = std::make_unique<VertexCompatibleWithBeam>(VertexDistanceXY(), algoconf.getParameter<double>("maxDistanceToBeam"));
    algorithm.is4D = algoconf.getParameter<bool>("is4D");

    // configure separate vertex time reconstruction if applicable
    // note that the vertex time could, in principle, also come from the clusterizer or the vertex fit
    if (algorithm.is4D) {
      auto const& pv_time_conf = algoconf.getParameter<edm::ParameterSet>("vertexTimeParameters");
      auto const vertexTimeAlgorithm = pv_time_conf.getParameter<std::string>("algorithm");
      LogDebug("PrimaryVertexProducer") << " vertexTimeParamers found  " << algorithm.label << " : [" << vertexTimeAlgorithm << "]";
      if (vertexTimeAlgorithm == "legacy4D") {
        algorithm.pv_time_estimator = std::make_unique<VertexTimeAlgorithmLegacy4D>(pv_time_conf.getParameter<edm::ParameterSet>("legacy4D"), consumesCollector());
        useTransientTrackTime = true;
      } else if (vertexTimeAlgorithm == "fromTracksPID") {
        algorithm.pv_time_estimator = std::make_unique<VertexTimeAlgorithmFromTracksPID>(pv_time_conf.getParameter<edm::ParameterSet>("fromTracksPID"), consumesCollector());
      } else {
        edm::LogWarning("PrimaryVertexProducer") << "unknown vertexTimeParameters.algorithm" << vertexTimeAlgorithm;
      }
    }

    produces<reco::VertexCollection>(algorithm.label);
  }

  if (useTransientTrackTime) {
    trkTimesToken = consumes(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
  }

  //check if this is a recovery iteration
  fRecoveryIteration = conf.getParameter<bool>("isRecoveryIteration");
  if (fRecoveryIteration) {
    if (algorithms.empty()) {
      throw VertexException("PrimaryVertexProducer: No algorithm specified.");
    } else if (algorithms.size() > 1) {
      throw VertexException("PrimaryVertexProducer: Running in Recovery mode and more than one algorithm specified. Please specify only one algorithm.");
    }
    recoveryVtxToken = consumes(conf.getParameter<edm::InputTag>("recoveryVtxCollection"));
  }
}

void PrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the BeamSpot, it will always be needed, even when not used as a constraint
  reco::BeamSpot beamSpot;
  bool validBS = true;

  auto const recoBeamSpotHandle = iEvent.getHandle(bsToken);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("PrimaryVertexProducer") << "Invalid BeamSpot: handle to reco::Beamspot product is not valid";
    validBS = false;
  }

  VertexState const beamVertexState{beamSpot};
  if (beamVertexState.error().cxx() <= 0. or beamVertexState.error().cyy() <= 0. or beamVertexState.error().czz() <= 0.) {
    edm::LogError("PrimaryVertexProducer") << "Invalid BeamSpot: non-positive beamspot errors " << beamVertexState.error().matrix();
    validBS = false;
  }

  // if this is a recovery iteration, check if we already have a valid PV
  if (fRecoveryIteration) {
    auto const& oldVertices = iEvent.get(recoveryVtxToken);
    // look for the first valid (not-BeamSpot) vertex
    for (auto const& old : oldVertices) {
      if (not old.isFake()) {
        // found a valid vertex, write the first one to the collection and return
        // otherwise continue with regular vertexing procedure
        auto result = std::make_unique<reco::VertexCollection>();
        result->emplace_back(old);
        iEvent.put(std::move(result), algorithms.begin()->label);
        return;
      }
    }
  }

  // initialise pv-time fitter
  for (auto& algo : algorithms) {
    if (algo.pv_time_estimator) {
      algo.pv_time_estimator->setEvent(iEvent, iSetup);
    }
  }

  // get RECO tracks from the event
  // tks can be used as a ptr to a reco::TrackCollection
  auto const tks = iEvent.getHandle(trkToken);

  // interface RECO tracks to vertex reconstruction
  const auto& theB = iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;

  if (useTransientTrackTime) {
    auto const& trackTimes = iEvent.get(trkTimesToken);
    auto const& trackTimeResos = iEvent.get(trkTimeResosToken);
    t_tks = theB.build(tks, beamSpot, trackTimes, trackTimeResos);
  } else {
    t_tks = theB.build(tks, beamSpot);
  }

  // select tracks
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);

  // clusterize tracks in Z
  std::vector<TransientVertex>&& clusters = theTrackClusterizer->vertices(seltks);

  if (fVerbose) {
    edm::LogPrint("PrimaryVertexProducer") << "Clustering returned " << clusters.size() << " clusters from " << seltks.size() << " selected tracks";
  }

  // vertex fits
  for (auto const& algo : algorithms) {

    std::vector<TransientVertex> pvs;
    if (algo.useBeamConstraint and not validBS){
      edm::LogError("PrimaryVertexProducer") << "Vertex Collection with label \"" << algo.label
          << "\" requires beam-constrained fit, but no valid BeamSpot in the Event. Returning empty collection of TransientVertexs.";
    } else {
      pvs = (algo.pv_fitter != nullptr) ? algo.pv_fitter->fit(seltks, clusters, beamSpot, algo.useBeamConstraint) : clusters;
    }

    // add vertex time
    if (algo.pv_time_estimator != nullptr) {
      algo.pv_time_estimator->fill_vertex_times(pvs);
    }

    // sort vertices by pt**2 vertex
    if (pvs.size() > 1) {
      std::sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }

    // output product (reco::VertexCollection)
    auto vColl = std::make_unique<reco::VertexCollection>();
    vColl->reserve(pvs.size());

    // select and convert transient vertices to (reco) vertices
    for (auto const& iv : pvs) {
      if(iv.isValid() and iv.degreesOfFreedom() >= algo.minNdof){
	reco::Vertex v = iv;
	if ((not validBS) or (*(algo.pv_selector))(v, beamVertexState)){
	  vColl->emplace_back(v);
	}
      }
    }

    if (fVerbose) {
      edm::LogPrint("PrimaryVertexProducer") << "PrimaryVertexProducer \"" << algo.label << "\" contains " << pvs.size() << " reco::Vertex candidates";
    }

    if (clusters.size() > 2 and clusters.size() > 2 * pvs.size()) {
      edm::LogWarning("PrimaryVertexProducer") << "More than 50% of candidate vertices lost (" << pvs.size() << " out of " << clusters.size() << ")";
    }

    if (pvs.empty() and seltks.size() > 5) {
      edm::LogWarning("PrimaryVertexProducer") << "No vertex found with " << seltks.size() << " tracks and " << clusters.size() << " vertex candidates";
    }

    if (vColl->empty()) {
      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl->emplace_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        edm::LogWarning("PrimaryVertexProducer") << "Zero recostructed vertices, will put reco::Vertex derived from dummy/fake BeamSpot into Event, BeamSpot has invalid errors: " << bse.matrix();
      } else {
        vColl->emplace_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        edm::LogWarning("PrimaryVertexProducer") << "Zero recostructed vertices, will put reco::Vertex derived from BeamSpot into Event.";
      }
    }

    if (fVerbose) {
      int ivtx = 0;
      for (auto const& vtx : *vColl) {
        edm::LogPrint("PrimaryVertexProducer") << "recvtx " << std::setw(3) << std::fixed << ivtx++ << " #trk " << std::setw(3) << vtx.tracksSize()
                  << " chi2 " << std::setw(5) << std::setprecision(1) << vtx.chi2() << " ndof " << std::setw(5)
                  << std::setprecision(1) << vtx.ndof() << " x " << std::setw(7) << std::setprecision(4)
                  << vtx.position().x() << " dx " << std::setw(6) << std::setprecision(4) << vtx.xError() << " y "
                  << std::setw(7) << std::setprecision(4) << vtx.position().y() << " dy " << std::setw(6)
                  << std::setprecision(4) << vtx.yError() << " z " << std::setw(8) << std::setprecision(4)
                 << vtx.position().z() << " dz " << std::setw(6) << std::setprecision(4) << vtx.zError()
 << " t " << std::setw(6) << std::setprecision(3) << vtx.t() << " dt " << std::setw(6)
                    << std::setprecision(3) << vtx.tError();
      }
    }

    iEvent.put(std::move(vColl), algo.label);
  }
}

void PrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription psd_pv_time;
  {
    edm::ParameterSetDescription psd1;
    VertexTimeAlgorithmLegacy4D::fillPSetDescription(psd1);
    psd_pv_time.add<edm::ParameterSetDescription>("legacy4D", psd1);

    edm::ParameterSetDescription psd2;
    VertexTimeAlgorithmFromTracksPID::fillPSetDescription(psd2);
    psd_pv_time.add<edm::ParameterSetDescription>("fromTracksPID", psd2);
  }
  psd_pv_time.add<std::string>("algorithm", "");  // default = none

  // vertex collections
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<std::string>("label", "");
    vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
    vpsd1.add<bool>("useBeamConstraint", false);
    vpsd1.add<double>("maxDistanceToBeam", 1.0);
    vpsd1.add<double>("chi2cutoff", 2.5);
    vpsd1.add<double>("mintrkweight", 0.0);
    vpsd1.add<double>("minNdof", 0.0);
    vpsd1.add<bool>("is4D", false);
    vpsd1.add<edm::ParameterSetDescription>("vertexTimeParameters", psd_pv_time);

    // two default values : with- and without beam constraint
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(2);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<std::string>("label", "");
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", false);
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("mintrkweight", 0.);
      temp2.addParameter<double>("minNdof", 0.0);
      temp2.addParameter<bool>("is4D", false);
      edm::ParameterSet temp_vertexTime;
      temp_vertexTime.addParameter<std::string>("algorithm", "");
      temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
      temp1.emplace_back(temp2);
    }
    {
      edm::ParameterSet temp2;
      temp2.addParameter<std::string>("label", "WithBS");
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", true);
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("mintrkweight", 0.);
      temp2.addParameter<double>("minNdof", 2.0);
      temp2.addParameter<bool>("is4D", false);
      edm::ParameterSet temp_vertexTime;
      temp_vertexTime.addParameter<std::string>("algorithm", "");
      temp2.addParameter<edm::ParameterSet>("vertexTimeParameters", temp_vertexTime);
      temp1.emplace_back(temp2);
    }
    desc.addVPSet("vertexCollections", vpsd1, temp1);
  }
  desc.addUntracked<bool>("verbose", false);
  {
    edm::ParameterSetDescription psd0;
    TrackFilterForPVFinding::fillPSetDescription(psd0);
    psd0.add<int>("numTracksThreshold", 0);  // HI only
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("TrackTimeResosLabel", edm::InputTag("dummy_default"));  // 4D only
  desc.add<edm::InputTag>("TrackTimesLabel", edm::InputTag("dummy_default"));      // 4D only

  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      DAClusterizerInZT_vect::fillPSetDescription(psd1);
      psd0.add<edm::ParameterSetDescription>("TkDAClusParameters", psd1);

      edm::ParameterSetDescription psd2;
      GapClusterizerInZ::fillPSetDescription(psd2);
      psd0.add<edm::ParameterSetDescription>("TkGapClusParameters", psd2);
    }
    psd0.add<std::string>("algorithm", "DA_vect");
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});
  desc.add<edm::ESInputTag>("transientTrackBuilder", edm::ESInputTag("", "TransientTrackBuilder"));

  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PrimaryVertexProducer);

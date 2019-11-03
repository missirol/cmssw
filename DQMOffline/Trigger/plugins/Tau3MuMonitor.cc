#include "DQMOffline/Trigger/plugins/Tau3MuMonitor.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

Tau3MuMonitor::Tau3MuMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("FolderName")),
      requireValidHLTPaths_(iConfig.getParameter<bool>("requireValidHLTPaths")),
      hltPathsAreValid_(false),
      tauToken_(mayConsume<reco::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
      pt_binning_(
          getHistoPSet(iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("ptPSet"))),
      eta_binning_(getHistoPSet(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("etaPSet"))),
      phi_binning_(getHistoPSet(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("phiPSet"))),
      mass_binning_(getHistoPSet(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet>("massPSet"))),
      genTriggerEventFlag_(new GenericTriggerEventFlag(
          iConfig.getParameter<edm::ParameterSet>("GenericTriggerEventPSet"), consumesCollector(), *this)) {
}

Tau3MuMonitor::~Tau3MuMonitor() throw() {

  if(genTriggerEventFlag_){ genTriggerEventFlag_.reset(); }
}

void Tau3MuMonitor::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const& iSetup) {

  // Initialize the GenericTriggerEventFlag
  if(genTriggerEventFlag_ && genTriggerEventFlag_->on()){ genTriggerEventFlag_->initRun(iRun, iSetup); }

  // check if every HLT path specified in numerator and denominator has a valid match in the HLT Menu
  hltPathsAreValid_ = (genTriggerEventFlag_ && genTriggerEventFlag_->on() && genTriggerEventFlag_->allHLTPathsAreValid());

  // if valid HLT paths are required,
  // create DQM outputs only if all paths are valid
  if (requireValidHLTPaths_ and (not hltPathsAreValid_)) {
    return;
  }

  std::string histname;

  std::string currentFolder = folderName_;
  ibooker.setCurrentFolder(currentFolder);

  // tau 3 mu 1D pt
  histname = "tau1DPt";
  tau1DPt_ = ibooker.book1D(histname, "", pt_binning_.nbins, pt_binning_.xmin, pt_binning_.xmax);
  tau1DPt_->setAxisTitle("3-#mu p_{T} [GeV]", 1);
  tau1DPt_->setAxisTitle("counts", 2);

  // tau 3 mu 1D eta
  histname = "tau1DEta";
  tau1DEta_ = ibooker.book1D(histname, "", eta_binning_.nbins, eta_binning_.xmin, eta_binning_.xmax);
  tau1DEta_->setAxisTitle("3-#mu #eta", 1);
  tau1DEta_->setAxisTitle("counts", 2);

  // tau 3 mu 1D phi
  histname = "tau1DPhi";
  tau1DPhi_ = ibooker.book1D(histname, "", phi_binning_.nbins, phi_binning_.xmin, phi_binning_.xmax);
  tau1DPhi_->setAxisTitle("3-#mu #phi", 1);
  tau1DPhi_->setAxisTitle("counts", 2);

  // tau 3 mu 1D mass
  histname = "tau1DMass";
  tau1DMass_ = ibooker.book1D(histname, "", mass_binning_.nbins, mass_binning_.xmin, mass_binning_.xmax);
  tau1DMass_->setAxisTitle("mass_{3#mu} [GeV]", 1);
  tau1DMass_->setAxisTitle("counts", 2);

  // tau 3 mu 2D eta vs phi
  histname = "tau2DEtaPhi";
  tau2DEtaPhi_ = ibooker.book2D(histname,
                                "",
                                eta_binning_.nbins,
                                eta_binning_.xmin,
                                eta_binning_.xmax,
                                phi_binning_.nbins,
                                phi_binning_.xmin,
                                phi_binning_.xmax);
  tau2DEtaPhi_->setAxisTitle("3-#mu #eta", 1);
  tau2DEtaPhi_->setAxisTitle("3-#mu #phi", 2);
}

void Tau3MuMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

  // if valid HLT paths are required,
  // analyze event only if all paths are valid
  if (requireValidHLTPaths_ and (not hltPathsAreValid_)) {
    return;
  }

  // require the trigger to be fired
  if (genTriggerEventFlag_->on() && !genTriggerEventFlag_->accept(iEvent, iSetup))
    return;

  // check if the previous event failed because of missing tau3mu collection.
  // Return silently, a warning must have been issued already at this point
  if (not validProduct_)
    return;

  // get ahold of the tau(3mu) collection
  edm::Handle<reco::CompositeCandidateCollection> tauHandle;
  iEvent.getByToken(tauToken_, tauHandle);

  // if the handle is not valid issue a warning (only for the forst occurrency)
  if (not tauHandle.isValid()) {

    edm::LogWarning("ProductNotValid") << "Tau3Mu trigger product not valid";
    validProduct_ = false;
    return;
  }

  // loop and fill
  for (auto const& itau : *tauHandle) {

    tau1DPt_->Fill(itau.pt());
    tau1DEta_->Fill(itau.eta());
    tau1DPhi_->Fill(itau.phi());
    tau1DMass_->Fill(itau.mass());
    tau2DEtaPhi_->Fill(itau.eta(), itau.phi());
  }
}

void Tau3MuMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("FolderName", "HLT/BPH/");
  desc.add<bool>("requireValidHLTPaths", false);

  desc.add<edm::InputTag>("taus", edm::InputTag("hltTauPt10MuPts511Mass1p2to2p3Iso", "Taus"));

  edm::ParameterSetDescription histoPSet;
  edm::ParameterSetDescription ptPSet;
  edm::ParameterSetDescription etaPSet;
  edm::ParameterSetDescription phiPSet;
  edm::ParameterSetDescription massPSet;
  fillHistoPSetDescription(ptPSet);
  fillHistoPSetDescription(etaPSet);
  fillHistoPSetDescription(phiPSet);
  fillHistoPSetDescription(massPSet);
  histoPSet.add<edm::ParameterSetDescription>("ptPSet", ptPSet);
  histoPSet.add<edm::ParameterSetDescription>("etaPSet", etaPSet);
  histoPSet.add<edm::ParameterSetDescription>("phiPSet", phiPSet);
  histoPSet.add<edm::ParameterSetDescription>("massPSet", massPSet);
  desc.add<edm::ParameterSetDescription>("histoPSet", histoPSet);

  edm::ParameterSetDescription genericTriggerEventPSet;
  genericTriggerEventPSet.add<bool>("andOr");
  genericTriggerEventPSet.add<edm::InputTag>("dcsInputTag", edm::InputTag("scalersRawToDigi"));
  genericTriggerEventPSet.add<std::vector<int> >("dcsPartitions", {});
  genericTriggerEventPSet.add<bool>("andOrDcs", false);
  genericTriggerEventPSet.add<bool>("errorReplyDcs", true);
  genericTriggerEventPSet.add<std::string>("dbLabel", "");
  genericTriggerEventPSet.add<bool>("andOrHlt", true);
  genericTriggerEventPSet.add<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults::HLT"));
  genericTriggerEventPSet.add<std::vector<std::string> >("hltPaths", {});
  genericTriggerEventPSet.add<std::string>("hltDBKey", "");
  genericTriggerEventPSet.add<bool>("errorReplyHlt", false);
  genericTriggerEventPSet.add<unsigned int>("verbosityLevel", 0);
  desc.add<edm::ParameterSetDescription>("GenericTriggerEventPSet", genericTriggerEventPSet);

  descriptions.add("tau3muMonitoring", desc);
}

DEFINE_FWK_MODULE(Tau3MuMonitor);

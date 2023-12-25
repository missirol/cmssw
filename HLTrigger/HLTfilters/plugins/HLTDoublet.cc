/** \class HLTDoublet
 *
 * See header file for documentation
 *
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "HLTDoublet.h"

#include <cmath>

template <typename T1, typename T2>
HLTDoublet<T1, T2>::HLTDoublet(const edm::ParameterSet& iConfig)
    : HLTFilter(iConfig),
      originTag1_(iConfig.template getParameter<std::vector<edm::InputTag>>("originTag1")),
      originTag2_(iConfig.template getParameter<std::vector<edm::InputTag>>("originTag2")),
      inputTag1_(iConfig.template getParameter<edm::InputTag>("inputTag1")),
      inputTag2_(iConfig.template getParameter<edm::InputTag>("inputTag2")),
      inputToken1_(consumes<trigger::TriggerFilterObjectWithRefs>(inputTag1_)),
      inputToken2_(consumes<trigger::TriggerFilterObjectWithRefs>(inputTag2_)),
      triggerType1_(iConfig.template getParameter<int>("triggerType1")),
      triggerType2_(iConfig.template getParameter<int>("triggerType2")),
      min_Dphi_(iConfig.template getParameter<double>("MinDphi")),
      max_Dphi_(iConfig.template getParameter<double>("MaxDphi")),
      min_Deta_(iConfig.template getParameter<double>("MinDeta")),
      max_Deta_(iConfig.template getParameter<double>("MaxDeta")),
      min_Minv_(iConfig.template getParameter<double>("MinMinv")),
      max_Minv_(iConfig.template getParameter<double>("MaxMinv")),
      min_DelR_(iConfig.template getParameter<double>("MinDelR")),
      max_DelR_(iConfig.template getParameter<double>("MaxDelR")),
      min_DelR2_(min_DelR_ * min_DelR_),
      max_DelR2_(max_DelR_ * max_DelR_),
      min_Pt_(iConfig.template getParameter<double>("MinPt")),
      max_Pt_(iConfig.template getParameter<double>("MaxPt")),
      min_N_(iConfig.template getParameter<int>("MinN")),
      same_(inputTag1_.encode() == inputTag2_.encode()),  // same collections to be compared?
      cutdphi_(min_Dphi_ <= max_Dphi_),                   // cut active?
      cutdeta_(min_Deta_ <= max_Deta_),                   // cut active?
      cutminv_(min_Minv_ <= max_Minv_),                   // cut active?
      cutdelr2_(min_DelR_ <= max_DelR_),                  // cut active?
      cutdelr2min_(min_DelR_ > 0),                        // cut active?
      cutdelr2max_(max_DelR_ > 0),                        // cut active?
      cutpt_(min_Pt_ <= max_Pt_)                          // cut active?
{
  LogDebug("") << "InputTags and cuts : " << inputTag1_.encode() << " " << inputTag2_.encode() << triggerType1_ << " "
               << triggerType2_ << " Dphi [" << min_Dphi_ << " " << max_Dphi_ << "]"
               << " Deta [" << min_Deta_ << " " << max_Deta_ << "]"
               << " Minv [" << min_Minv_ << " " << max_Minv_ << "]"
               << " DelR [" << min_DelR_ << " " << max_DelR_ << "]"
               << " Pt   [" << min_Pt_ << " " << max_Pt_ << "]"
               << " MinN=" << min_N_ << " same=" << same_ << " dphi=" << cutdphi_ << " deta=" << cutdeta_
               << " minv=" << cutminv_ << " delr=" << cutdelr2_ << " delrmin=" << cutdelr2min_
               << " delrmax=" << cutdelr2max_ << " pt=" << cutpt_;
}

template <typename T1, typename T2>
HLTDoublet<T1, T2>::~HLTDoublet() = default;
template <typename T1, typename T2>
void HLTDoublet<T1, T2>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  std::vector<edm::InputTag> originTag1(1, edm::InputTag("hltOriginal1"));
  std::vector<edm::InputTag> originTag2(1, edm::InputTag("hltOriginal2"));
  desc.add<std::vector<edm::InputTag>>("originTag1", originTag1);
  desc.add<std::vector<edm::InputTag>>("originTag2", originTag2);
  desc.add<edm::InputTag>("inputTag1", edm::InputTag("hltFiltered1"));
  desc.add<edm::InputTag>("inputTag2", edm::InputTag("hltFiltered22"));
  desc.add<int>("triggerType1", 0);
  desc.add<int>("triggerType2", 0);
  desc.add<double>("MinDphi", +1.0);
  desc.add<double>("MaxDphi", -1.0);
  desc.add<double>("MinDeta", +1.0);
  desc.add<double>("MaxDeta", -1.0);
  desc.add<double>("MinMinv", +1.0);
  desc.add<double>("MaxMinv", -1.0);
  desc.add<double>("MinDelR", +1.0);
  desc.add<double>("MaxDelR", -1.0);
  desc.add<double>("MinPt", +1.0);
  desc.add<double>("MaxPt", -1.0);
  desc.add<int>("MinN", 1);
  descriptions.add(defaultModuleLabel<HLTDoublet<T1, T2>>(), desc);
}

//
// member functions
//

// ------------ method called to produce the data  ------------
template <typename T1, typename T2>
bool HLTDoublet<T1, T2>::hltFilter(edm::Event& iEvent,
                                   const edm::EventSetup& iSetup,
                                   trigger::TriggerFilterObjectWithRefs& filterproduct) const {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  // All HLT filters must create and fill an HLT filter object,
  // recording any reconstructed physics objects satisfying (or not)
  // this HLT filter, and place it in the Event.

  bool accept(false);

  LogVerbatim("HLTDoublet") << " XXX " << moduleLabel() << " 0 ";

  std::vector<T1Ref> coll1;
  std::vector<T2Ref> coll2;

  // get hold of pre-filtered object collections
  Handle<TriggerFilterObjectWithRefs> handle1, handle2;
  if (iEvent.getByToken(inputToken1_, handle1) && iEvent.getByToken(inputToken2_, handle2)) {
    handle1->getObjects(triggerType1_, coll1);
    handle2->getObjects(triggerType2_, coll2);
    const size_type n1(coll1.size());
    const size_type n2(coll2.size());

    if (saveTags()) {
      InputTag tagOld;
      for (unsigned int i = 0; i < originTag1_.size(); ++i) {
        filterproduct.addCollectionTag(originTag1_[i]);
        LogVerbatim("HLTDoublet") << " XXX " << moduleLabel() << " 1a/" << i << " " << originTag1_[i].encode();
      }
      tagOld = InputTag();
      for (size_type i1 = 0; i1 != n1; ++i1) {
        const ProductID pid(coll1[i1].id());
        const auto& prov = iEvent.getStableProvenance(pid);
        const string& label(prov.moduleLabel());
        const string& instance(prov.productInstanceName());
        const string& process(prov.processName());
        InputTag tagNew(InputTag(label, instance, process));
        if (tagOld.encode() != tagNew.encode()) {
          filterproduct.addCollectionTag(tagNew);
          tagOld = tagNew;
          LogVerbatim("HLTDoublet") << " XXX " << moduleLabel() << " 1b " << tagNew.encode();
        }
      }
      for (unsigned int i = 0; i < originTag2_.size(); ++i) {
        filterproduct.addCollectionTag(originTag2_[i]);
        LogVerbatim("HLTDoublet") << " XXX " << moduleLabel() << " 2a/" << i << " " << originTag2_[i].encode();
      }
      tagOld = InputTag();
      for (size_type i2 = 0; i2 != n2; ++i2) {
        const ProductID pid(coll2[i2].id());
        const auto& prov = iEvent.getStableProvenance(pid);
        const string& label(prov.moduleLabel());
        const string& instance(prov.productInstanceName());
        const string& process(prov.processName());
        InputTag tagNew(InputTag(label, instance, process));
        if (tagOld.encode() != tagNew.encode()) {
          filterproduct.addCollectionTag(tagNew);
          tagOld = tagNew;
          LogVerbatim("HLTDoublet") << " XXX " << moduleLabel() << " 2b " << tagNew.encode();
        }
      }
    }

    int n(0);
    T1Ref r1;
    T2Ref r2;
    Particle::LorentzVector p1, p2, p;
    for (unsigned int i1 = 0; i1 != n1; i1++) {
      r1 = coll1[i1];
      p1 = r1->p4();
      unsigned int I(0);
      if (same_) {
        I = i1 + 1;
      }
      for (unsigned int i2 = I; i2 != n2; i2++) {
        r2 = coll2[i2];
        p2 = r2->p4();

        double const Deta(std::abs(p1.eta() - p2.eta()));
        if (cutdeta_ and (min_Deta_ > Deta or Deta > max_Deta_))
          continue;

        double const Dphi = std::abs(reco::deltaPhi(p1.phi(), p2.phi()));
        if (cutdphi_ and (min_Dphi_ > Dphi or Dphi > max_Dphi_))
          continue;

        double const DelR2 = Deta * Deta + Dphi * Dphi;
        if (cutdelr2_ and cutdelr2min_ and min_DelR2_ > DelR2)
          continue;
        if (cutdelr2_ and cutdelr2max_ and max_DelR2_ < DelR2)
          continue;

        p = p1 + p2;

        double const Pt = p.pt();
        if (cutpt_ and (min_Pt_ > Pt or Pt > max_Pt_))
          continue;

        double const Minv = std::abs(p.mass());
        if (cutminv_ and (min_Minv_ > Minv or Minv > max_Minv_))
          continue;

        n++;
        filterproduct.addObject(triggerType1_, r1);
        filterproduct.addObject(triggerType2_, r2);
      }
    }

    // filter decision
    accept = (n >= min_N_);
  }

  return accept;
}

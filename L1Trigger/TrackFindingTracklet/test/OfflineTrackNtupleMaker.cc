// ---------------------------------------------------------------------------
// OfflineTrackNtupleMaker
//
// Offline tracks (generalTracks) matched to TrackingParticles, in a flat tree
// for measuring offline tracking efficiency, fake rate and duplicate rate.
//
//   trk_*       offline tracks     trk_isTrue, trk_matchtp_*, trk_matchtp_quality
//   tp_*        truth particles    tp_nmatch
//   matchtrk_*  best track per TP, index-aligned with the tp_* vectors
//
// The three quantities come straight out of these:
//
//   efficiency     tp_nmatch > 0        over all tp   (differentially in pt/eta/phi)
//   fake rate      trk_isTrue == 0      over all trk
//   duplicate rate tp_nmatch > 1        over tp with tp_nmatch > 0
//
// Branch names deliberately mirror the L1 ntuple (trk_/tp_/matchtrk_) so the
// analysis macros in analyses/tracking_performance port over with little change.
//
// Matching is hit-based and done by the standard chain, which this plugin only
// consumes -- it redoes no matching itself:
//   tpClusterProducer -> quickTrackAssociatorByHits
//                     -> trackingParticleRecoTrackAsssociation
//
// This is an offline-only job: it needs no L1 re-emulation and runs directly on
// step3.root. The L1 ntuplizers are untouched.
// ---------------------------------------------------------------------------

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include <TTree.h>
#include <cmath>
#include <vector>

class OfflineTrackNtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit OfflineTrackNtupleMaker(const edm::ParameterSet&);
  ~OfflineTrackNtupleMaker() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ---- configuration ----
  const edm::EDGetTokenT<edm::View<reco::Track>> trackToken_;
  const edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
  const edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimToken_;
  const edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const double tpMinPt_;
  const double tpMaxEta_;
  const double tpMaxZ0_;
  const int tpMinNHits_;
  const bool onlyChargedTP_;

  TTree* tree_ = nullptr;

  // ---- reco track branches ----
  std::vector<float> *m_trk_pt = nullptr, *m_trk_eta = nullptr, *m_trk_phi = nullptr;
  std::vector<float> *m_trk_d0 = nullptr, *m_trk_z0 = nullptr;
  std::vector<float> *m_trk_chi2 = nullptr, *m_trk_chi2dof = nullptr;
  std::vector<int> *m_trk_nhit = nullptr, *m_trk_nlayer = nullptr;
  std::vector<int> *m_trk_npixhit = nullptr, *m_trk_charge = nullptr;
  std::vector<int> *m_trk_algo = nullptr, *m_trk_highPurity = nullptr;
  // truth for this track: quality is the association purity (fraction of shared hits)
  std::vector<int> *m_trk_isTrue = nullptr;
  std::vector<float> *m_trk_matchtp_pt = nullptr, *m_trk_matchtp_eta = nullptr;
  std::vector<float> *m_trk_matchtp_phi = nullptr, *m_trk_matchtp_z0 = nullptr;
  std::vector<int> *m_trk_matchtp_pdgid = nullptr;
  std::vector<float> *m_trk_matchtp_quality = nullptr;

  // ---- tracking particle branches ----
  std::vector<float> *m_tp_pt = nullptr, *m_tp_eta = nullptr, *m_tp_phi = nullptr;
  std::vector<float> *m_tp_z0 = nullptr, *m_tp_d0 = nullptr, *m_tp_lxy = nullptr;
  std::vector<int> *m_tp_pdgid = nullptr, *m_tp_charge = nullptr;
  std::vector<int> *m_tp_nhit = nullptr, *m_tp_ngenpart = nullptr;
  // number of offline tracks matched to this TP (0 = missed, >1 = duplicates)
  std::vector<int> *m_tp_nmatch = nullptr;

  // ---- best matched tracks per TP, index-aligned with the tp_ vectors ----
  std::vector<float> *m_matchtrk_pt = nullptr, *m_matchtrk_eta = nullptr;
  std::vector<float> *m_matchtrk_phi = nullptr;
  std::vector<float> *m_matchtrk_d0 = nullptr, *m_matchtrk_z0 = nullptr;
  std::vector<float> *m_matchtrk_chi2dof = nullptr, *m_matchtrk_quality = nullptr;
  std::vector<int> *m_matchtrk_nhit = nullptr;

  std::vector<std::vector<float>*> allF_;
  std::vector<std::vector<int>*> allI_;

  void bookF(const char* name, std::vector<float>** v) {
    *v = new std::vector<float>;
    allF_.push_back(*v);
    tree_->Branch(name, *v);
  }
  void bookI(const char* name, std::vector<int>** v) {
    *v = new std::vector<int>;
    allI_.push_back(*v);
    tree_->Branch(name, *v);
  }
};

OfflineTrackNtupleMaker::OfflineTrackNtupleMaker(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<edm::View<reco::Track>>(iConfig.getParameter<edm::InputTag>("TrackInputTag"))),
      tpToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag"))),
      recoToSimToken_(consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("AssociatorInputTag"))),
      simToRecoToken_(consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("AssociatorInputTag"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotInputTag"))),
      tpMinPt_(iConfig.getParameter<double>("TP_minPt")),
      tpMaxEta_(iConfig.getParameter<double>("TP_maxEta")),
      tpMaxZ0_(iConfig.getParameter<double>("TP_maxZ0")),
      tpMinNHits_(iConfig.getParameter<int>("TP_minNHits")),
      onlyChargedTP_(iConfig.getParameter<bool>("TP_onlyCharged")) {
  usesResource("TFileService");
}

void OfflineTrackNtupleMaker::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("eventTree", "offline track / tracking particle matching");

  bookF("trk_pt", &m_trk_pt);
  bookF("trk_eta", &m_trk_eta);
  bookF("trk_phi", &m_trk_phi);
  bookF("trk_d0", &m_trk_d0);
  bookF("trk_z0", &m_trk_z0);
  bookF("trk_chi2", &m_trk_chi2);
  bookF("trk_chi2dof", &m_trk_chi2dof);
  bookI("trk_nhit", &m_trk_nhit);
  bookI("trk_nlayer", &m_trk_nlayer);
  bookI("trk_npixhit", &m_trk_npixhit);
  bookI("trk_charge", &m_trk_charge);
  bookI("trk_algo", &m_trk_algo);
  bookI("trk_highPurity", &m_trk_highPurity);
  bookI("trk_isTrue", &m_trk_isTrue);
  bookF("trk_matchtp_pt", &m_trk_matchtp_pt);
  bookF("trk_matchtp_eta", &m_trk_matchtp_eta);
  bookF("trk_matchtp_phi", &m_trk_matchtp_phi);
  bookF("trk_matchtp_z0", &m_trk_matchtp_z0);
  bookI("trk_matchtp_pdgid", &m_trk_matchtp_pdgid);
  bookF("trk_matchtp_quality", &m_trk_matchtp_quality);

  bookF("tp_pt", &m_tp_pt);
  bookF("tp_eta", &m_tp_eta);
  bookF("tp_phi", &m_tp_phi);
  bookF("tp_z0", &m_tp_z0);
  bookF("tp_d0", &m_tp_d0);
  bookF("tp_lxy", &m_tp_lxy);
  bookI("tp_pdgid", &m_tp_pdgid);
  bookI("tp_charge", &m_tp_charge);
  bookI("tp_nmatch", &m_tp_nmatch);
  bookI("tp_nhit", &m_tp_nhit);
  bookI("tp_ngenpart", &m_tp_ngenpart);

  bookF("matchtrk_pt", &m_matchtrk_pt);
  bookF("matchtrk_eta", &m_matchtrk_eta);
  bookF("matchtrk_phi", &m_matchtrk_phi);
  bookF("matchtrk_d0", &m_matchtrk_d0);
  bookF("matchtrk_z0", &m_matchtrk_z0);
  bookF("matchtrk_chi2dof", &m_matchtrk_chi2dof);
  bookF("matchtrk_quality", &m_matchtrk_quality);
  bookI("matchtrk_nhit", &m_matchtrk_nhit);
}

void OfflineTrackNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  for (auto* v : allF_) v->clear();
  for (auto* v : allI_) v->clear();

  edm::Handle<edm::View<reco::Track>> tracks;
  iEvent.getByToken(trackToken_, tracks);
  edm::Handle<TrackingParticleCollection> tps;
  iEvent.getByToken(tpToken_, tps);
  edm::Handle<reco::RecoToSimCollection> recoToSim;
  iEvent.getByToken(recoToSimToken_, recoToSim);
  edm::Handle<reco::SimToRecoCollection> simToReco;
  iEvent.getByToken(simToRecoToken_, simToReco);
  edm::Handle<reco::BeamSpot> bs;
  iEvent.getByToken(beamSpotToken_, bs);

  const math::XYZPoint bsPos =
      bs.isValid() ? math::XYZPoint(bs->x0(), bs->y0(), bs->z0()) : math::XYZPoint(0, 0, 0);

  // Products may legitimately be absent (e.g. a file without the association
  // chain run over it). Fill what is there rather than throwing; empty trk_*
  // branches then point straight at a missing association.
  const bool haveOffline = tracks.isValid() && recoToSim.isValid() && simToReco.isValid();
  const bool haveTP = tps.isValid();

  // ---------------------------------------------------------------- reco -> sim
  for (size_t i = 0; haveOffline && i < tracks->size(); ++i) {
    edm::RefToBase<reco::Track> trkRef(tracks, i);
    const reco::Track& trk = *trkRef;

    m_trk_pt->push_back(trk.pt());
    m_trk_eta->push_back(trk.eta());
    m_trk_phi->push_back(trk.phi());
    m_trk_d0->push_back(trk.dxy(bsPos));
    m_trk_z0->push_back(trk.dz(bsPos));
    m_trk_chi2->push_back(trk.chi2());
    m_trk_chi2dof->push_back(trk.ndof() > 0 ? trk.chi2() / trk.ndof() : -999.);
    m_trk_nhit->push_back(trk.numberOfValidHits());
    m_trk_nlayer->push_back(trk.hitPattern().trackerLayersWithMeasurement());
    m_trk_npixhit->push_back(trk.hitPattern().numberOfValidPixelHits());
    m_trk_charge->push_back(trk.charge());
    m_trk_algo->push_back(trk.algo());
    m_trk_highPurity->push_back(trk.quality(reco::TrackBase::highPurity) ? 1 : 0);

    // The associator returns the matched TPs ordered by quality; take the best.
    int isTrue = 0;
    float mpt = -999, meta = -999, mphi = -999, mz0 = -999, mq = -999;
    int mpdg = 0;
    auto found = recoToSim->find(trkRef);
    if (found != recoToSim->end() && !found->val.empty()) {
      const auto& best = found->val.front();
      const TrackingParticleRef& tp = best.first;
      isTrue = 1;
      mpt = tp->pt();
      meta = tp->eta();
      mphi = tp->phi();
      mz0 = tp->vz();
      mpdg = tp->pdgId();
      mq = best.second;  // association quality / purity
    }
    m_trk_isTrue->push_back(isTrue);
    m_trk_matchtp_pt->push_back(mpt);
    m_trk_matchtp_eta->push_back(meta);
    m_trk_matchtp_phi->push_back(mphi);
    m_trk_matchtp_z0->push_back(mz0);
    m_trk_matchtp_pdgid->push_back(mpdg);
    m_trk_matchtp_quality->push_back(mq);
  }

  // ---------------------------------------------------------------- sim -> reco
  for (size_t i = 0; haveTP && i < tps->size(); ++i) {
    TrackingParticleRef tpRef(tps, i);  // the associator works in Refs
    const TrackingParticle& tp = *tpRef;

    // Selection kept deliberately loose and configurable; TP_minNHits = 0 and
    // TP_minPt = 0 give the full stored collection.
    if (onlyChargedTP_ && tp.charge() == 0) continue;
    if (tp.pt() < tpMinPt_) continue;
    if (std::abs(tp.eta()) > tpMaxEta_) continue;
    if (std::abs(tp.vz()) > tpMaxZ0_) continue;
    if (tpMinNHits_ > 0 && tp.numberOfTrackerHits() < tpMinNHits_) continue;

    const float lxy = std::sqrt(tp.vx() * tp.vx() + tp.vy() * tp.vy());
    // d0 of the production vertex about the beamspot, sign convention as in the
    // L1 ntuplizer.
    const float d0 = -(tp.vx() - bsPos.x()) * std::sin(tp.phi()) +
                     (tp.vy() - bsPos.y()) * std::cos(tp.phi());

    m_tp_pt->push_back(tp.pt());
    m_tp_eta->push_back(tp.eta());
    m_tp_phi->push_back(tp.phi());
    m_tp_z0->push_back(tp.vz());
    m_tp_d0->push_back(d0);
    m_tp_lxy->push_back(lxy);
    m_tp_pdgid->push_back(tp.pdgId());
    m_tp_charge->push_back(tp.charge());
    m_tp_nhit->push_back(tp.numberOfTrackerHits());
    // 0 => produced by GEANT rather than by the event generator
    m_tp_ngenpart->push_back(static_cast<int>(tp.genParticles().size()));

    // --- offline match for this TP ---
    int nmatch = 0;
    float bpt = -999, beta = -999, bphi = -999, bd0 = -999, bz0 = -999;
    float bchi2dof = -999, bq = -999;
    int bnhit = -1;
    if (haveOffline) {
      auto found = simToReco->find(tpRef);
      if (found != simToReco->end() && !found->val.empty()) {
      nmatch = found->val.size();
      const auto& best = found->val.front();  // highest quality first
      const reco::Track& t = *best.first;
      bpt = t.pt();
      beta = t.eta();
      bphi = t.phi();
      bd0 = t.dxy(bsPos);
      bz0 = t.dz(bsPos);
      bchi2dof = t.ndof() > 0 ? t.chi2() / t.ndof() : -999.;
      bnhit = t.numberOfValidHits();
        bq = best.second;
      }
    }
    m_tp_nmatch->push_back(nmatch);
    m_matchtrk_pt->push_back(bpt);
    m_matchtrk_eta->push_back(beta);
    m_matchtrk_phi->push_back(bphi);
    m_matchtrk_d0->push_back(bd0);
    m_matchtrk_z0->push_back(bz0);
    m_matchtrk_chi2dof->push_back(bchi2dof);
    m_matchtrk_nhit->push_back(bnhit);
    m_matchtrk_quality->push_back(bq);
  }

  tree_->Fill();
}

void OfflineTrackNtupleMaker::endJob() {
  edm::LogVerbatim("OfflineTrackNtupleMaker") << "OfflineTrackNtupleMaker: done";
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(OfflineTrackNtupleMaker);

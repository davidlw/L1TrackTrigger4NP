//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

////////////////////
// Added by Wei Li for hits of modules
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTStubAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "L1Trigger/TrackFindingTracklet/interface/HitPatternHelper.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CLHEP/Units/PhysicalConstants.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;

class L1TrackHitNtupleMaker : public one::EDAnalyzer<one::WatchRuns, one::SharedResources> {
public:
  explicit L1TrackHitNtupleMaker(const edm::ParameterSet& iConfig);
  ~L1TrackHitNtupleMaker() override;

  void beginJob() override;
  void endJob() override;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void beginRun(const Run& iEvent, const EventSetup& iSetup) override {}
  void endRun(const Run& iEvent, const EventSetup& iSetup) override {}

protected:
private:
  edm::ParameterSet config;

  int MyProcess;       
  bool DebugMode;      
  bool SaveAllTracks;  
  bool SaveStubs;      
  int L1Tk_nPar;       
  int TP_minNStub;  
  int TP_minNStubLayer;  
  double TP_minPt;       
  double TP_maxEta;      
  double TP_maxZ0;       
  int L1Tk_minNStub;     

  bool TrackingInJets;  

  edm::InputTag L1TrackInputTag;       
  edm::InputTag MCTruthTrackInputTag;  
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag L1ClusterInputTag;     
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  edm::InputTag TrackingVertexInputTag;
  edm::InputTag GenJetInputTag;

  edm::EDGetTokenT<edmNew::DetSetVector<Phase2TrackerCluster1D>> phase2OTClustersToken_; 

  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> ttClusterToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>> ttStubToken_;
  edm::EDGetTokenT<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> ttClusterMCTruthToken_;
  edm::EDGetTokenT<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> ttStubMCTruthToken_;

  edm::EDGetTokenT<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> ttTrackToken_;
  edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> ttTrackMCTruthToken_;

  edm::EDGetTokenT<std::vector<TrackingParticle>> TrackingParticleToken_;
  edm::EDGetTokenT<std::vector<TrackingVertex>> TrackingVertexToken_;

  edm::EDGetTokenT<std::vector<reco::GenJet>> GenJetToken_;

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> getTokenTrackerGeom_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> getTokenTrackerTopo_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> getTokenBField_;
  edm::ESGetToken<hph::Setup, hph::SetupRcd> getTokenHPHSetup_;
  edm::ESGetToken<tt::Setup, tt::SetupRcd> getTokenSetup_;
  edm::ESGetToken<trackerTFP::LayerEncoding, trackerTFP::DataFormatsRcd> getTokenLayerEncoding_;

  bool available_;  

  TTree* eventTree;

  // all L1 tracks
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<float>* m_trk_d0;  
  std::vector<float>* m_trk_z0;
  std::vector<float>* m_trk_chi2;
  std::vector<float>* m_trk_chi2_dof;
  std::vector<float>* m_trk_chi2rphi;
  std::vector<float>* m_trk_chi2rphi_dof;
  std::vector<float>* m_trk_chi2rz;
  std::vector<float>* m_trk_chi2rz_dof;
  std::vector<float>* m_trk_bendchi2;
  std::vector<int>* m_trk_nstub;
  std::vector<int>* m_trk_lhits;
  std::vector<int>* m_trk_dhits;
  std::vector<int>* m_trk_seed;
  std::vector<int>* m_trk_hitpattern;
  std::vector<int>* m_trk_lhits_hitpattern;  
  std::vector<int>* m_trk_dhits_hitpattern;  
  std::vector<int>* m_trk_nPSstub_hitpattern;
  std::vector<int>* m_trk_n2Sstub_hitpattern;
  std::vector<int>* m_trk_nLostPSstub_hitpattern;
  std::vector<int>* m_trk_nLost2Sstub_hitpattern;
  std::vector<int>* m_trk_nLoststub_V1_hitpattern;  
  std::vector<int>* m_trk_nLoststub_V2_hitpattern;  
  std::vector<int>* m_trk_charge;
  std::vector<unsigned int>* m_trk_phiSector;
  std::vector<int>* m_trk_etaSector;
  std::vector<int>* m_trk_genuine;
  std::vector<int>* m_trk_loose;
  std::vector<int>* m_trk_unknown;
  std::vector<int>* m_trk_combinatoric;
  std::vector<int>* m_trk_fake;  
  std::vector<float>* m_trk_MVA1;
  std::vector<int>* m_trk_matchtp_pdgid;
  std::vector<float>* m_trk_matchtp_pt;
  std::vector<float>* m_trk_matchtp_eta;
  std::vector<float>* m_trk_matchtp_phi;
  std::vector<float>* m_trk_matchtp_z0;
  std::vector<float>* m_trk_matchtp_lxy;
  std::vector<float>* m_trk_matchtp_d0;
  std::vector<int>* m_trk_injet;          
  std::vector<int>* m_trk_injet_highpt;   
  std::vector<int>* m_trk_injet_vhighpt;  
  std::vector<std::vector<int>>* m_trk_layers;

  // all tracking particles
  std::vector<float>* m_tp_pt;
  std::vector<float>* m_tp_eta;
  std::vector<float>* m_tp_phi;
  std::vector<float>* m_tp_lxy;
  std::vector<float>* m_tp_d0;
  std::vector<float>* m_tp_z0;
  std::vector<float>* m_tp_d0_prod;
  std::vector<float>* m_tp_z0_prod;
  std::vector<int>* m_tp_pdgid;
  std::vector<int>* m_tp_nmatch;
  std::vector<int>* m_tp_nstub;
  std::vector<int>* m_tp_eventid;
  std::vector<int>* m_tp_charge;
  std::vector<int>* m_tp_injet;
  std::vector<int>* m_tp_injet_highpt;
  std::vector<int>* m_tp_injet_vhighpt;

  // *L1 track* properties if m_tp_nmatch > 0
  std::vector<float>* m_matchtrk_pt;
  std::vector<float>* m_matchtrk_eta;
  std::vector<float>* m_matchtrk_phi;
  std::vector<float>* m_matchtrk_d0;  
  std::vector<float>* m_matchtrk_z0;
  std::vector<float>* m_matchtrk_chi2;
  std::vector<float>* m_matchtrk_chi2_dof;
  std::vector<float>* m_matchtrk_chi2rphi;
  std::vector<float>* m_matchtrk_chi2rphi_dof;
  std::vector<float>* m_matchtrk_chi2rz;
  std::vector<float>* m_matchtrk_chi2rz_dof;
  std::vector<float>* m_matchtrk_bendchi2;
  std::vector<float>* m_matchtrk_MVA1;
  std::vector<int>* m_matchtrk_nstub;
  std::vector<int>* m_matchtrk_lhits;
  std::vector<int>* m_matchtrk_dhits;
  std::vector<int>* m_matchtrk_seed;
  std::vector<int>* m_matchtrk_hitpattern;
  std::vector<int>* m_matchtrk_charge;
  std::vector<int>* m_matchtrk_injet;
  std::vector<int>* m_matchtrk_injet_highpt;
  std::vector<int>* m_matchtrk_injet_vhighpt;

  // ALL stubs
  std::vector<float>* m_allstub_x;
  std::vector<float>* m_allstub_y;
  std::vector<float>* m_allstub_z;

  std::vector<int>* m_allstub_isBarrel;  
  std::vector<int>* m_allstub_layer;
  std::vector<int>* m_allstub_isPSmodule;
  std::vector<int>* m_allstub_isTiltedBarrel;

  std::vector<float>* m_allstub_trigDisplace;
  std::vector<float>* m_allstub_trigOffset;
  std::vector<float>* m_allstub_trigPos;
  std::vector<float>* m_allstub_trigBend;

  // CLUSTER BRANCHES (Wei Li / Rice)
  std::vector<float> *cluster_x, *cluster_y, *cluster_z; // Added cluster_z
  std::vector<int> *cluster_layer, *cluster_isBarrel, *cluster_halfModule, *cluster_isPS, *cluster_chipId, *cluster_sensor;
  std::vector<uint32_t> *cluster_detid;

  // INCLUSIVE TTCLUSTER BRANCHES
  std::vector<float>* m_ttclus_x;
  std::vector<float>* m_ttclus_y;
  std::vector<float>* m_ttclus_z;
  std::vector<int>* m_ttclus_layer;
  std::vector<int>* m_ttclus_isBarrel;
  std::vector<int>* m_ttclus_isPS;
  std::vector<int>* m_ttclus_width;
  std::vector<int>* m_ttclus_sensor; 

  // stub associated with tracking particle ?
  std::vector<int>* m_allstub_matchTP_pdgid;  
  std::vector<float>* m_allstub_matchTP_pt;   
  std::vector<float>* m_allstub_matchTP_eta;  
  std::vector<float>* m_allstub_matchTP_phi;  

  std::vector<int>* m_allstub_genuine;

  // track jet variables
  std::vector<float>* m_jet_eta;
  std::vector<float>* m_jet_phi;
  std::vector<float>* m_jet_pt;
  std::vector<float>* m_jet_tp_sumpt;
  std::vector<float>* m_jet_trk_sumpt;
  std::vector<float>* m_jet_matchtrk_sumpt;
};

//////////////
// CONSTRUCTOR
L1TrackHitNtupleMaker::L1TrackHitNtupleMaker(edm::ParameterSet const& iConfig) : config(iConfig) {
  usesResource("TFileService");
  MyProcess = iConfig.getParameter<int>("MyProcess");
  DebugMode = iConfig.getParameter<bool>("DebugMode");
  SaveAllTracks = iConfig.getParameter<bool>("SaveAllTracks");
  SaveStubs = iConfig.getParameter<bool>("SaveStubs");
  L1Tk_nPar = iConfig.getParameter<int>("L1Tk_nPar");
  TP_minNStub = iConfig.getParameter<int>("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter<int>("TP_minNStubLayer");
  TP_minPt = iConfig.getParameter<double>("TP_minPt");
  TP_maxEta = iConfig.getParameter<double>("TP_maxEta");
  TP_maxZ0 = iConfig.getParameter<double>("TP_maxZ0");
  L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub = iConfig.getParameter<int>("L1Tk_minNStub");

  TrackingInJets = iConfig.getParameter<bool>("TrackingInJets");

  L1StubInputTag = iConfig.getParameter<edm::InputTag>("L1StubInputTag");
  L1ClusterInputTag = iConfig.getParameter<edm::InputTag>("L1ClusterInputTag"); 
  MCTruthClusterInputTag = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  MCTruthStubInputTag = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  TrackingVertexInputTag = iConfig.getParameter<edm::InputTag>("TrackingVertexInputTag");
  GenJetInputTag = iConfig.getParameter<edm::InputTag>("GenJetInputTag");

  phase2OTClustersToken_ = consumes<edmNew::DetSetVector<Phase2TrackerCluster1D>>(iConfig.getParameter<edm::InputTag>("phase2OTClusters")); 

  ttTrackToken_ = consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>>(MCTruthTrackInputTag);
  ttStubToken_ = consumes<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>>(L1StubInputTag);
  ttClusterToken_ = consumes<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>>(L1ClusterInputTag);
  ttClusterMCTruthToken_ = consumes<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>>(MCTruthClusterInputTag);
  ttStubMCTruthToken_ = consumes<TTStubAssociationMap<Ref_Phase2TrackerDigi_>>(MCTruthStubInputTag);

  TrackingParticleToken_ = consumes<std::vector<TrackingParticle>>(TrackingParticleInputTag);
  TrackingVertexToken_ = consumes<std::vector<TrackingVertex>>(TrackingVertexInputTag);
  GenJetToken_ = consumes<std::vector<reco::GenJet>>(GenJetInputTag);

  getTokenTrackerGeom_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  getTokenTrackerTopo_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
  getTokenBField_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  getTokenHPHSetup_ = esConsumes<hph::Setup, hph::SetupRcd>();
  getTokenSetup_ = esConsumes<tt::Setup, tt::SetupRcd>();
  getTokenLayerEncoding_ = esConsumes<trackerTFP::LayerEncoding, trackerTFP::DataFormatsRcd>();
}

/////////////
// DESTRUCTOR
L1TrackHitNtupleMaker::~L1TrackHitNtupleMaker() {}

//////////
// END JOB
void L1TrackHitNtupleMaker::endJob() {
  edm::LogVerbatim("Tracklet") << "L1TrackHitNtupleMaker::endJob";

  delete m_trk_pt;
  delete m_trk_eta;
  delete m_trk_phi;
  delete m_trk_z0;
  delete m_trk_d0;
  delete m_trk_chi2;
  delete m_trk_chi2_dof;
  delete m_trk_chi2rphi;
  delete m_trk_chi2rphi_dof;
  delete m_trk_chi2rz;
  delete m_trk_chi2rz_dof;
  delete m_trk_bendchi2;
  delete m_trk_nstub;
  delete m_trk_lhits;
  delete m_trk_dhits;
  delete m_trk_seed;
  delete m_trk_hitpattern;
  delete m_trk_lhits_hitpattern;
  delete m_trk_dhits_hitpattern;
  delete m_trk_nPSstub_hitpattern;
  delete m_trk_n2Sstub_hitpattern;
  delete m_trk_nLostPSstub_hitpattern;
  delete m_trk_nLost2Sstub_hitpattern;
  delete m_trk_nLoststub_V1_hitpattern;
  delete m_trk_nLoststub_V2_hitpattern;
  delete m_trk_charge;
  delete m_trk_phiSector;
  delete m_trk_etaSector;
  delete m_trk_genuine;
  delete m_trk_loose;
  delete m_trk_unknown;
  delete m_trk_combinatoric;
  delete m_trk_fake;
  delete m_trk_MVA1;
  delete m_trk_matchtp_pdgid;
  delete m_trk_matchtp_pt;
  delete m_trk_matchtp_eta;
  delete m_trk_matchtp_phi;
  delete m_trk_matchtp_z0;
  delete m_trk_matchtp_lxy;
  delete m_trk_matchtp_d0;
  delete m_trk_injet;
  delete m_trk_injet_highpt;
  delete m_trk_injet_vhighpt;
  delete m_trk_layers;

  delete m_tp_pt;
  delete m_tp_eta;
  delete m_tp_phi;
  delete m_tp_lxy;
  delete m_tp_d0;
  delete m_tp_z0;
  delete m_tp_d0_prod;
  delete m_tp_z0_prod;
  delete m_tp_pdgid;
  delete m_tp_nmatch;
  delete m_tp_nstub;
  delete m_tp_eventid;
  delete m_tp_charge;
  delete m_tp_injet;
  delete m_tp_injet_highpt;
  delete m_tp_injet_vhighpt;

  delete m_matchtrk_pt;
  delete m_matchtrk_eta;
  delete m_matchtrk_phi;
  delete m_matchtrk_z0;
  delete m_matchtrk_d0;
  delete m_matchtrk_chi2;
  delete m_matchtrk_chi2_dof;
  delete m_matchtrk_chi2rphi;
  delete m_matchtrk_chi2rphi_dof;
  delete m_matchtrk_chi2rz;
  delete m_matchtrk_chi2rz_dof;
  delete m_matchtrk_bendchi2;
  delete m_matchtrk_MVA1;
  delete m_matchtrk_nstub;
  delete m_matchtrk_dhits;
  delete m_matchtrk_lhits;
  delete m_matchtrk_seed;
  delete m_matchtrk_hitpattern;
  delete m_matchtrk_charge;
  delete m_matchtrk_injet;
  delete m_matchtrk_injet_highpt;
  delete m_matchtrk_injet_vhighpt;

  delete m_allstub_x;
  delete m_allstub_y;
  delete m_allstub_z;
  delete m_allstub_isBarrel;
  delete m_allstub_layer;
  delete m_allstub_isPSmodule;
  delete m_allstub_isTiltedBarrel;
  delete m_allstub_trigDisplace;
  delete m_allstub_trigOffset;
  delete m_allstub_trigPos;
  delete m_allstub_trigBend;
  delete m_allstub_matchTP_pdgid;
  delete m_allstub_matchTP_pt;
  delete m_allstub_matchTP_eta;
  delete m_allstub_matchTP_phi;
  delete m_allstub_genuine;

  delete cluster_x;
  delete cluster_y;
  delete cluster_z; // Clean up cluster_z
  delete cluster_layer;
  delete cluster_isBarrel;
  delete cluster_halfModule;
  delete cluster_isPS;
  delete cluster_chipId;
  delete cluster_sensor;
  delete cluster_detid;

  delete m_ttclus_x;
  delete m_ttclus_y;
  delete m_ttclus_z;
  delete m_ttclus_layer;
  delete m_ttclus_isBarrel;
  delete m_ttclus_isPS;
  delete m_ttclus_width;
  delete m_ttclus_sensor;

  delete m_jet_eta;
  delete m_jet_phi;
  delete m_jet_pt;
  delete m_jet_tp_sumpt;
  delete m_jet_trk_sumpt;
  delete m_jet_matchtrk_sumpt;
}

////////////
// BEGIN JOB
void L1TrackHitNtupleMaker::beginJob() {
  edm::LogVerbatim("Tracklet") << "L1TrackHitNtupleMaker::beginJob";

  edm::Service<TFileService> fs;
  available_ = fs.isAvailable();
  if (not available_)
    return;  

  m_trk_pt = new std::vector<float>;
  m_trk_eta = new std::vector<float>;
  m_trk_phi = new std::vector<float>;
  m_trk_z0 = new std::vector<float>;
  m_trk_d0 = new std::vector<float>;
  m_trk_chi2 = new std::vector<float>;
  m_trk_chi2_dof = new std::vector<float>;
  m_trk_chi2rphi = new std::vector<float>;
  m_trk_chi2rphi_dof = new std::vector<float>;
  m_trk_chi2rz = new std::vector<float>;
  m_trk_chi2rz_dof = new std::vector<float>;
  m_trk_bendchi2 = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_lhits = new std::vector<int>;
  m_trk_dhits = new std::vector<int>;
  m_trk_seed = new std::vector<int>;
  m_trk_hitpattern = new std::vector<int>;
  m_trk_lhits_hitpattern = new std::vector<int>;
  m_trk_dhits_hitpattern = new std::vector<int>;
  m_trk_nPSstub_hitpattern = new std::vector<int>;
  m_trk_n2Sstub_hitpattern = new std::vector<int>;
  m_trk_nLostPSstub_hitpattern = new std::vector<int>;
  m_trk_nLost2Sstub_hitpattern = new std::vector<int>;
  m_trk_nLoststub_V1_hitpattern = new std::vector<int>;
  m_trk_nLoststub_V2_hitpattern = new std::vector<int>;
  m_trk_charge = new std::vector<int>;
  m_trk_phiSector = new std::vector<unsigned int>;
  m_trk_etaSector = new std::vector<int>;
  m_trk_genuine = new std::vector<int>;
  m_trk_loose = new std::vector<int>;
  m_trk_unknown = new std::vector<int>;
  m_trk_combinatoric = new std::vector<int>;
  m_trk_fake = new std::vector<int>;
  m_trk_MVA1 = new std::vector<float>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt = new std::vector<float>;
  m_trk_matchtp_eta = new std::vector<float>;
  m_trk_matchtp_phi = new std::vector<float>;
  m_trk_matchtp_z0 = new std::vector<float>;
  m_trk_matchtp_lxy = new std::vector<float>;
  m_trk_matchtp_d0 = new std::vector<float>;
  m_trk_injet = new std::vector<int>;
  m_trk_injet_highpt = new std::vector<int>;
  m_trk_injet_vhighpt = new std::vector<int>;
  m_trk_layers = new std::vector<std::vector<int>>;

  m_tp_pt = new std::vector<float>;
  m_tp_eta = new std::vector<float>;
  m_tp_phi = new std::vector<float>;
  m_tp_lxy = new std::vector<float>;
  m_tp_d0 = new std::vector<float>;
  m_tp_z0 = new std::vector<float>;
  m_tp_d0_prod = new std::vector<float>;
  m_tp_z0_prod = new std::vector<float>;
  m_tp_pdgid = new std::vector<int>;
  m_tp_nmatch = new std::vector<int>;
  m_tp_nstub = new std::vector<int>;
  m_tp_eventid = new std::vector<int>;
  m_tp_charge = new std::vector<int>;
  m_tp_injet = new std::vector<int>;
  m_tp_injet_highpt = new std::vector<int>;
  m_tp_injet_vhighpt = new std::vector<int>;

  m_matchtrk_pt = new std::vector<float>;
  m_matchtrk_eta = new std::vector<float>;
  m_matchtrk_phi = new std::vector<float>;
  m_matchtrk_z0 = new std::vector<float>;
  m_matchtrk_d0 = new std::vector<float>;
  m_matchtrk_chi2 = new std::vector<float>;
  m_matchtrk_chi2_dof = new std::vector<float>;
  m_matchtrk_chi2rphi = new std::vector<float>;
  m_matchtrk_chi2rphi_dof = new std::vector<float>;
  m_matchtrk_chi2rz = new std::vector<float>;
  m_matchtrk_chi2rz_dof = new std::vector<float>;
  m_matchtrk_bendchi2 = new std::vector<float>;
  m_matchtrk_MVA1 = new std::vector<float>;
  m_matchtrk_nstub = new std::vector<int>;
  m_matchtrk_dhits = new std::vector<int>;
  m_matchtrk_lhits = new std::vector<int>;
  m_matchtrk_seed = new std::vector<int>;
  m_matchtrk_hitpattern = new std::vector<int>;
  m_matchtrk_charge = new std::vector<int>;
  m_matchtrk_injet = new std::vector<int>;
  m_matchtrk_injet_highpt = new std::vector<int>;
  m_matchtrk_injet_vhighpt = new std::vector<int>;

  m_allstub_x = new std::vector<float>;
  m_allstub_y = new std::vector<float>;
  m_allstub_z = new std::vector<float>;
  m_allstub_isBarrel = new std::vector<int>;
  m_allstub_layer = new std::vector<int>;
  m_allstub_isPSmodule = new std::vector<int>;
  m_allstub_isTiltedBarrel = new std::vector<int>;
  m_allstub_trigDisplace = new std::vector<float>;
  m_allstub_trigOffset = new std::vector<float>;
  m_allstub_trigPos = new std::vector<float>;
  m_allstub_trigBend = new std::vector<float>;
  m_allstub_matchTP_pdgid = new std::vector<int>;
  m_allstub_matchTP_pt = new std::vector<float>;
  m_allstub_matchTP_eta = new std::vector<float>;
  m_allstub_matchTP_phi = new std::vector<float>;
  m_allstub_genuine = new std::vector<int>;

  cluster_x = new std::vector<float>; cluster_y = new std::vector<float>; cluster_z = new std::vector<float>; // Initialize cluster_z
  cluster_layer = new std::vector<int>; cluster_isBarrel = new std::vector<int>;
  cluster_halfModule = new std::vector<int>; cluster_detid = new std::vector<uint32_t>;
  cluster_isPS = new std::vector<int>; cluster_chipId = new std::vector<int>; cluster_sensor = new std::vector<int>;

  m_ttclus_x        = new std::vector<float>;
  m_ttclus_y        = new std::vector<float>;
  m_ttclus_z        = new std::vector<float>;
  m_ttclus_layer    = new std::vector<int>;
  m_ttclus_isBarrel = new std::vector<int>;
  m_ttclus_isPS     = new std::vector<int>;
  m_ttclus_width    = new std::vector<int>;
  m_ttclus_sensor   = new std::vector<int>;

  m_jet_eta = new std::vector<float>;
  m_jet_phi = new std::vector<float>;
  m_jet_pt = new std::vector<float>;
  m_jet_tp_sumpt = new std::vector<float>;
  m_jet_trk_sumpt = new std::vector<float>;
  m_jet_matchtrk_sumpt = new std::vector<float>;

  eventTree = fs->make<TTree>("eventTree", "Event tree");

  if (SaveAllTracks) {
    eventTree->Branch("trk_pt", &m_trk_pt);
    eventTree->Branch("trk_eta", &m_trk_eta);
    eventTree->Branch("trk_phi", &m_trk_phi);
    eventTree->Branch("trk_d0", &m_trk_d0);
    eventTree->Branch("trk_z0", &m_trk_z0);
    eventTree->Branch("trk_chi2", &m_trk_chi2);
    eventTree->Branch("trk_chi2_dof", &m_trk_chi2_dof);
    eventTree->Branch("trk_chi2rphi", &m_trk_chi2rphi);
    eventTree->Branch("trk_chi2rphi_dof", &m_trk_chi2rphi_dof);
    eventTree->Branch("trk_chi2rz", &m_trk_chi2rz);
    eventTree->Branch("trk_chi2rz_dof", &m_trk_chi2rz_dof);
    eventTree->Branch("trk_bendchi2", &m_trk_bendchi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);
    eventTree->Branch("trk_lhits", &m_trk_lhits);
    eventTree->Branch("trk_dhits", &m_trk_dhits);
    eventTree->Branch("trk_seed", &m_trk_seed);
    eventTree->Branch("trk_hitpattern", &m_trk_hitpattern);
    eventTree->Branch("trk_lhits_hitpattern", &m_trk_lhits_hitpattern);
    eventTree->Branch("trk_dhits_hitpattern", &m_trk_dhits_hitpattern);
    eventTree->Branch("trk_nPSstub_hitpattern", &m_trk_nPSstub_hitpattern);
    eventTree->Branch("trk_n2Sstub_hitpattern", &m_trk_n2Sstub_hitpattern);
    eventTree->Branch("trk_nLostPSstub_hitpattern", &m_trk_nLostPSstub_hitpattern);
    eventTree->Branch("trk_n2Sstub_hitpattern", &m_trk_n2Sstub_hitpattern);
    eventTree->Branch("trk_nLost2Sstub_hitpattern", &m_trk_nLost2Sstub_hitpattern);
    eventTree->Branch("trk_nLoststub_V1_hitpattern", &m_trk_nLoststub_V1_hitpattern);
    eventTree->Branch("trk_nLoststub_V2_hitpattern", &m_trk_nLoststub_V2_hitpattern);
    eventTree->Branch("trk_charge", &m_trk_charge);
    eventTree->Branch("trk_phiSector", &m_trk_phiSector);
    eventTree->Branch("trk_etaSector", &m_trk_etaSector);
    eventTree->Branch("trk_genuine", &m_trk_genuine);
    eventTree->Branch("trk_loose", &m_trk_loose);
    eventTree->Branch("trk_unknown", &m_trk_unknown);
    eventTree->Branch("trk_combinatoric", &m_trk_combinatoric);
    eventTree->Branch("trk_fake", &m_trk_fake);
    eventTree->Branch("trk_MVA1", &m_trk_MVA1);
    eventTree->Branch("trk_matchtp_pdgid", &m_trk_matchtp_pdgid);
    eventTree->Branch("trk_matchtp_pt", &m_trk_matchtp_pt);
    eventTree->Branch("trk_matchtp_eta", &m_trk_matchtp_eta);
    eventTree->Branch("trk_matchtp_phi", &m_trk_matchtp_phi);
    eventTree->Branch("trk_matchtp_z0", &m_trk_matchtp_z0);
    eventTree->Branch("trk_matchtp_lxy", &m_trk_matchtp_lxy);
    eventTree->Branch("trk_matchtp_d0", &m_trk_matchtp_d0);
    if (TrackingInJets) {
      eventTree->Branch("trk_injet", &m_trk_injet);
      eventTree->Branch("trk_injet_highpt", &m_trk_injet_highpt);
      eventTree->Branch("trk_injet_vhighpt", &m_trk_injet_vhighpt);
    }
    eventTree->Branch("m_trk_layers", &m_trk_layers);
  }

  eventTree->Branch("tp_pt", &m_tp_pt);
  eventTree->Branch("tp_eta", &m_tp_eta);
  eventTree->Branch("tp_phi", &m_tp_phi);
  eventTree->Branch("tp_lxy", &m_tp_lxy);
  eventTree->Branch("tp_d0", &m_tp_d0);
  eventTree->Branch("tp_z0", &m_tp_z0);
  eventTree->Branch("tp_d0_prod", &m_tp_d0_prod);
  eventTree->Branch("tp_z0_prod", &m_tp_z0_prod);
  eventTree->Branch("tp_pdgid", &m_tp_pdgid);
  eventTree->Branch("tp_nmatch", &m_tp_nmatch);
  eventTree->Branch("tp_nstub", &m_tp_nstub);
  eventTree->Branch("tp_eventid", &m_tp_eventid);
  eventTree->Branch("tp_charge", &m_tp_charge);
  if (TrackingInJets) {
    eventTree->Branch("tp_injet", &m_tp_injet);
    eventTree->Branch("tp_injet_highpt", &m_tp_injet_highpt);
    eventTree->Branch("tp_injet_vhighpt", &m_tp_injet_vhighpt);
  }

  eventTree->Branch("matchtrk_pt", &m_matchtrk_pt);
  eventTree->Branch("matchtrk_eta", &m_matchtrk_eta);
  eventTree->Branch("matchtrk_phi", &m_matchtrk_phi);
  eventTree->Branch("matchtrk_z0", &m_matchtrk_z0);
  eventTree->Branch("matchtrk_d0", &m_matchtrk_d0);
  eventTree->Branch("matchtrk_chi2", &m_matchtrk_chi2);
  eventTree->Branch("matchtrk_chi2_dof", &m_matchtrk_chi2_dof);
  eventTree->Branch("matchtrk_chi2rphi", &m_matchtrk_chi2rphi);
  eventTree->Branch("matchtrk_chi2rphi_dof", &m_matchtrk_chi2rphi_dof);
  eventTree->Branch("matchtrk_chi2rz", &m_matchtrk_chi2rz);
  eventTree->Branch("matchtrk_chi2rz_dof", &m_matchtrk_chi2rz_dof);
  eventTree->Branch("matchtrk_bendchi2", &m_matchtrk_bendchi2);
  eventTree->Branch("matchtrk_MVA1", &m_matchtrk_MVA1);
  eventTree->Branch("matchtrk_nstub", &m_matchtrk_nstub);
  eventTree->Branch("matchtrk_lhits", &m_matchtrk_lhits);
  eventTree->Branch("matchtrk_dhits", &m_matchtrk_dhits);
  eventTree->Branch("matchtrk_seed", &m_matchtrk_seed);
  eventTree->Branch("matchtrk_hitpattern", &m_matchtrk_hitpattern);
  eventTree->Branch("matchtrk_charge", &m_matchtrk_charge);
  if (TrackingInJets) {
    eventTree->Branch("matchtrk_injet", &m_matchtrk_injet);
    eventTree->Branch("matchtrk_injet_highpt", &m_matchtrk_injet_highpt);
    eventTree->Branch("matchtrk_injet_vhighpt", &m_matchtrk_injet_vhighpt);
  }

  if (SaveStubs) {
    eventTree->Branch("allstub_x", &m_allstub_x);
    eventTree->Branch("allstub_y", &m_allstub_y);
    eventTree->Branch("allstub_z", &m_allstub_z);
    eventTree->Branch("allstub_isBarrel", &m_allstub_isBarrel);
    eventTree->Branch("allstub_layer", &m_allstub_layer);
    eventTree->Branch("allstub_isPSmodule", &m_allstub_isPSmodule);
    eventTree->Branch("allstub_isTiltedBarrel", &m_allstub_isTiltedBarrel);
    eventTree->Branch("allstub_trigDisplace", &m_allstub_trigDisplace);
    eventTree->Branch("allstub_trigOffset", &m_allstub_trigOffset);
    eventTree->Branch("allstub_trigPos", &m_allstub_trigPos);
    eventTree->Branch("allstub_trigBend", &m_allstub_trigBend);
    eventTree->Branch("allstub_matchTP_pdgid", &m_allstub_matchTP_pdgid);
    eventTree->Branch("allstub_matchTP_pt", &m_allstub_matchTP_pt);
    eventTree->Branch("allstub_matchTP_eta", &m_allstub_matchTP_eta);
    eventTree->Branch("allstub_matchTP_phi", &m_allstub_matchTP_phi);
    eventTree->Branch("allstub_genuine", &m_allstub_genuine);

    // Book Wei offline branches (with cluster_z included)
    eventTree->Branch("cluster_x", &cluster_x); eventTree->Branch("cluster_y", &cluster_y); eventTree->Branch("cluster_z", &cluster_z);
    eventTree->Branch("cluster_layer", &cluster_layer); eventTree->Branch("cluster_isBarrel", &cluster_isBarrel);
    eventTree->Branch("cluster_halfModule", &cluster_halfModule); eventTree->Branch("cluster_detid", &cluster_detid);
    eventTree->Branch("cluster_isPS", &cluster_isPS); eventTree->Branch("cluster_chipId", &cluster_chipId);
    eventTree->Branch("cluster_sensor", &cluster_sensor);

    eventTree->Branch("ttclus_x", &m_ttclus_x);
    eventTree->Branch("ttclus_y", &m_ttclus_y);
    eventTree->Branch("ttclus_z", &m_ttclus_z);
    eventTree->Branch("ttclus_layer", &m_ttclus_layer);
    eventTree->Branch("ttclus_isBarrel", &m_ttclus_isBarrel);
    eventTree->Branch("ttclus_isPS", &m_ttclus_isPS);
    eventTree->Branch("ttclus_width", &m_ttclus_width);
    eventTree->Branch("ttclus_sensor", &m_ttclus_sensor);
  }

  if (TrackingInJets) {
    eventTree->Branch("jet_eta", &m_jet_eta);
    eventTree->Branch("jet_phi", &m_jet_phi);
    eventTree->Branch("jet_pt", &m_jet_pt);
    eventTree->Branch("jet_tp_sumpt", &m_jet_tp_sumpt);
    eventTree->Branch("jet_trk_sumpt", &m_jet_trk_sumpt);
    eventTree->Branch("jet_matchtrk_sumpt", &m_jet_matchtrk_sumpt);
  }
}

//////////
// ANALYZE
void L1TrackHitNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (not available_)
    return;  

  if (!(MyProcess == 13 || MyProcess == 11 || MyProcess == 211 || MyProcess == 6 || MyProcess == 15 || MyProcess == 1)) {
    edm::LogVerbatim("Tracklet") << "The specified MyProcess is invalid! Exiting...";
    return;
  }

  if (!(L1Tk_nPar == 4 || L1Tk_nPar == 5)) {
    edm::LogVerbatim("Tracklet") << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting...";
    return;
  }

  if (SaveAllTracks) {
    m_trk_pt->clear(); m_trk_eta->clear(); m_trk_phi->clear(); m_trk_d0->clear(); m_trk_z0->clear();
    m_trk_chi2->clear(); m_trk_chi2_dof->clear(); m_trk_chi2rphi->clear(); m_trk_chi2rphi_dof->clear();
    m_trk_chi2rz->clear(); m_trk_chi2rz_dof->clear(); m_trk_bendchi2->clear(); m_trk_nstub->clear();
    m_trk_lhits->clear(); m_trk_dhits->clear(); m_trk_seed->clear(); m_trk_hitpattern->clear();
    m_trk_lhits_hitpattern->clear(); m_trk_dhits_hitpattern->clear(); m_trk_nPSstub_hitpattern->clear();
    m_trk_n2Sstub_hitpattern->clear(); m_trk_nLostPSstub_hitpattern->clear(); m_trk_nLost2Sstub_hitpattern->clear();
    m_trk_nLoststub_V1_hitpattern->clear(); m_trk_nLoststub_V2_hitpattern->clear(); m_trk_charge->clear();
    m_trk_phiSector->clear(); m_trk_etaSector->clear(); m_trk_genuine->clear(); m_trk_loose->clear();
    m_trk_unknown->clear(); m_trk_combinatoric->clear(); m_trk_fake->clear(); m_trk_MVA1->clear();
    m_trk_matchtp_pdgid->clear(); m_trk_matchtp_pt->clear(); m_trk_matchtp_eta->clear();
    m_trk_matchtp_phi->clear(); m_trk_matchtp_z0->clear(); m_trk_matchtp_lxy->clear();
    m_trk_matchtp_d0->clear(); m_trk_injet->clear(); m_trk_injet_highpt->clear();
    m_trk_injet_vhighpt->clear(); m_trk_layers->clear();
  }

  m_tp_pt->clear(); m_tp_eta->clear(); m_tp_phi->clear(); m_tp_lxy->clear(); m_tp_d0->clear(); m_tp_z0->clear();
  m_tp_d0_prod->clear(); m_tp_z0_prod->clear(); m_tp_pdgid->clear(); m_tp_nmatch->clear(); m_tp_nstub->clear();
  m_tp_eventid->clear(); m_tp_charge->clear(); m_tp_injet->clear(); m_tp_injet_highpt->clear(); m_tp_injet_vhighpt->clear();

  m_matchtrk_pt->clear(); m_matchtrk_eta->clear(); m_matchtrk_phi->clear(); m_matchtrk_z0->clear();
  m_matchtrk_d0->clear(); m_matchtrk_chi2->clear(); m_matchtrk_chi2_dof->clear(); m_matchtrk_chi2rphi->clear();
  m_matchtrk_chi2rphi_dof->clear(); m_matchtrk_chi2rz->clear(); m_matchtrk_chi2rz_dof->clear();
  m_matchtrk_bendchi2->clear(); m_matchtrk_MVA1->clear(); m_matchtrk_nstub->clear(); m_matchtrk_lhits->clear();
  m_matchtrk_dhits->clear(); m_matchtrk_seed->clear(); m_matchtrk_hitpattern->clear(); m_matchtrk_charge->clear();
  m_matchtrk_injet->clear(); m_matchtrk_injet_highpt->clear(); m_matchtrk_injet_vhighpt->clear();

  if (SaveStubs) {
    m_allstub_x->clear(); m_allstub_y->clear(); m_allstub_z->clear();
    m_allstub_isBarrel->clear(); m_allstub_layer->clear(); m_allstub_isPSmodule->clear(); m_allstub_isTiltedBarrel->clear();
    m_allstub_trigDisplace->clear(); m_allstub_trigOffset->clear(); m_allstub_trigPos->clear(); m_allstub_trigBend->clear();
    m_allstub_matchTP_pdgid->clear(); m_allstub_matchTP_pt->clear(); m_allstub_matchTP_eta->clear(); m_allstub_matchTP_phi->clear();
    m_allstub_genuine->clear();

    cluster_x->clear(); cluster_y->clear(); cluster_z->clear(); // Clear cluster_z
    cluster_layer->clear(); cluster_isBarrel->clear();
    cluster_halfModule->clear(); cluster_detid->clear(); cluster_isPS->clear(); cluster_chipId->clear();
    cluster_sensor->clear();

    m_ttclus_x->clear(); m_ttclus_y->clear(); m_ttclus_z->clear();
    m_ttclus_layer->clear(); m_ttclus_isBarrel->clear(); m_ttclus_isPS->clear();
    m_ttclus_width->clear(); m_ttclus_sensor->clear();
  }

  m_jet_eta->clear(); m_jet_phi->clear(); m_jet_pt->clear();
  m_jet_tp_sumpt->clear(); m_jet_trk_sumpt->clear(); m_jet_matchtrk_sumpt->clear();

  edm::Handle<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  edm::Handle<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>> TTStubHandle;
  if (SaveStubs)
    iEvent.getByToken(ttStubToken_, TTStubHandle);

  edm::Handle<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
  edm::Handle<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  edm::Handle<std::vector<TrackingParticle>> TrackingParticleHandle;
  edm::Handle<std::vector<TrackingVertex>> TrackingVertexHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle = iSetup.getHandle(getTokenTrackerGeom_);
  edm::ESHandle<TrackerTopology> tTopoHandle = iSetup.getHandle(getTokenTrackerTopo_);
  edm::ESHandle<MagneticField> bFieldHandle = iSetup.getHandle(getTokenBField_);
  edm::ESHandle<hph::Setup> hphHandle = iSetup.getHandle(getTokenHPHSetup_);
  edm::ESHandle<tt::Setup> handleSetup = iSetup.getHandle(getTokenSetup_);
  edm::ESHandle<trackerTFP::LayerEncoding> handleLayerEncoding = iSetup.getHandle(getTokenLayerEncoding_);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();
  const hph::Setup* hphSetup = hphHandle.product();
  const tt::Setup* setup = handleSetup.product();
  const trackerTFP::LayerEncoding* layerEncoding = handleLayerEncoding.product();

  // ----------------------------------------------------------------------------------------------
  // loop over L1 stubs
  // ----------------------------------------------------------------------------------------------
  if (SaveStubs) {
    for (auto gd = theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {
      DetId detid = (*gd)->geographicalId();
      if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
        continue;
      if (!tTopo->isLower(detid))
        continue;              
      DetId stackDetid = tTopo->stack(detid);  

      if (TTStubHandle->find(stackDetid) == TTStubHandle->end())
        continue;

      edmNew::DetSet<TTStub<Ref_Phase2TrackerDigi_>> stubs = (*TTStubHandle)[stackDetid];
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit(detid);
      const auto* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
      const PixelTopology* topol = dynamic_cast<const PixelTopology*>(&(theGeomDet->specificTopology()));

      for (auto stubIter = stubs.begin(); stubIter != stubs.end(); ++stubIter) {
        edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>> tempStubPtr =
            edmNew::makeRefTo(TTStubHandle, stubIter);

        int isBarrel = 0; int layer = -1;
        if (detid.subdetId() == StripSubdetector::TOB) {
          isBarrel = 1; layer = static_cast<int>(tTopo->layer(detid));
        } else if (detid.subdetId() == StripSubdetector::TID) {
          isBarrel = 0; layer = static_cast<int>(tTopo->layer(detid));
        }

        int isPSmodule = (topol->nrows() == 960) ? 1 : 0;
        const unsigned int tobSide = tTopo->tobSide(detid);  
        int isTiltedBarrel = (isBarrel == 1 && (tobSide == 1 || tobSide == 2)) ? 1 : 0;

        MeasurementPoint coords = tempStubPtr->clusterRef(0)->findAverageLocalCoordinatesCentered();
        LocalPoint clustlp = topol->localPosition(coords);
        GlobalPoint posStub = theGeomDet->surface().toGlobal(clustlp);

        m_allstub_x->push_back(posStub.x()); m_allstub_y->push_back(posStub.y()); m_allstub_z->push_back(posStub.z());
        m_allstub_isBarrel->push_back(isBarrel); m_allstub_layer->push_back(layer);
        m_allstub_isPSmodule->push_back(isPSmodule); m_allstub_isTiltedBarrel->push_back(isTiltedBarrel);
        m_allstub_trigDisplace->push_back(tempStubPtr->rawBend()); m_allstub_trigOffset->push_back(tempStubPtr->bendOffset());
        m_allstub_trigPos->push_back(tempStubPtr->innerClusterPosition()); m_allstub_trigBend->push_back(tempStubPtr->bendFE());

        edm::Ptr<TrackingParticle> my_tp = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);
        int myTP_pdgid = -999; float myTP_pt = -999; float myTP_eta = -999; float myTP_phi = -999;

        if (!my_tp.isNull()) {
          if (my_tp->eventId().event() == 0) {
            myTP_pdgid = my_tp->pdgId(); myTP_pt = my_tp->p4().pt();
            myTP_eta = my_tp->p4().eta(); myTP_phi = my_tp->p4().phi();
          }
        }
        m_allstub_matchTP_pdgid->push_back(myTP_pdgid); m_allstub_matchTP_pt->push_back(myTP_pt);
        m_allstub_matchTP_eta->push_back(myTP_eta); m_allstub_matchTP_phi->push_back(myTP_phi);
        m_allstub_genuine->push_back(MCTruthTTStubHandle->isGenuine(tempStubPtr) ? 1 : 0);
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  // Wei Li Offline Clusters
  // -----------------------------------------------------------------------------------------------
  edm::Handle<edmNew::DetSetVector<Phase2TrackerCluster1D>> phase2OTClusters;
  iEvent.getByToken(phase2OTClustersToken_, phase2OTClusters);

  if (phase2OTClusters.isValid()) {
    for (auto const& detSet : *phase2OTClusters) {
      DetId detId(detSet.detId());
      if (detId.det() != DetId::Tracker) continue;

      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit(detId);
      const auto* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
      if (!theGeomDet) continue;
      const PixelTopology* topol = dynamic_cast<const PixelTopology*>(&(theGeomDet->specificTopology()));

      int isPS = (topol->nrows() == 960) ? 1 : 0;
      int layer = tTopo->getOTLayerNumber(detId);
      int isBarrel = (detId.subdetId() == StripSubdetector::TOB);
      int isUpper = tTopo->isUpper(detId);

      for (auto const& cluster : detSet) {
        float center = cluster.center();
        int chipIdx = isPS ? (!isUpper ? static_cast<int>(center)/118 : static_cast<int>(center)/120) : static_cast<int>(center)/127;
        int halfMod = (center < 1024.0) ? 0 : 1;

        MeasurementPoint mp(center, 0.5);
        GlobalPoint gp = theGeomDet->surface().toGlobal(theGeomDet->topology().localPosition(mp));

        cluster_x->push_back(gp.x()); cluster_y->push_back(gp.y()); cluster_z->push_back(gp.z()); // Added cluster_z extraction
        cluster_layer->push_back(layer); cluster_isBarrel->push_back(isBarrel ? 1 : 0);
        cluster_halfModule->push_back(halfMod); cluster_detid->push_back(detId.rawId());
        cluster_isPS->push_back(isPS); cluster_chipId->push_back(chipIdx);
        cluster_sensor->push_back(isUpper ? 1 : 0);
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  // Inclusive Online TTClusters Extraction Implementation
  // -----------------------------------------------------------------------------------------------
  if (SaveStubs) {
    edm::Handle<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> TTClusterHandle;
    iEvent.getByToken(ttClusterToken_, TTClusterHandle);

    if (TTClusterHandle.isValid()) {
      for (auto const& detSet : *TTClusterHandle) {
        DetId detId(detSet.detId());
        if (detId.det() != DetId::Tracker) continue;

        const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit(detId);
        const auto* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
        if (!theGeomDet) continue;
        const PixelTopology* topol = dynamic_cast<const PixelTopology*>(&(theGeomDet->specificTopology()));

        int isPS = (topol->nrows() == 960) ? 1 : 0;
        int layer = tTopo->getOTLayerNumber(detId);
        int isBarrel = (detId.subdetId() == StripSubdetector::TOB) ? 1 : 0;

        for (auto clusIter = detSet.begin(); clusIter != detSet.end(); ++clusIter) {
          MeasurementPoint coords = clusIter->findAverageLocalCoordinatesCentered();
          GlobalPoint posClust = theGeomDet->surface().toGlobal(topol->localPosition(coords));

          m_ttclus_x->push_back(posClust.x());
          m_ttclus_y->push_back(posClust.y());
          m_ttclus_z->push_back(posClust.z());
          m_ttclus_layer->push_back(layer);
          m_ttclus_isBarrel->push_back(isBarrel);
          m_ttclus_isPS->push_back(isPS);
          m_ttclus_width->push_back(clusIter->getRows().size());
          m_ttclus_sensor->push_back(clusIter->getStackMember()); 
        }
      }
    }
  }

  // ----------------------------------------------------------------------------------------------
  // tracking in jets
  // ----------------------------------------------------------------------------------------------
  std::vector<math::XYZTLorentzVector> v_jets;
  std::vector<int> v_jets_highpt; std::vector<int> v_jets_vhighpt;

  if (TrackingInJets) {
    edm::Handle<std::vector<reco::GenJet>> GenJetHandle;
    iEvent.getByToken(GenJetToken_, GenJetHandle);

    if (GenJetHandle.isValid()) {
      for (auto iterGenJet = GenJetHandle->begin(); iterGenJet != GenJetHandle->end(); ++iterGenJet) {
        if (iterGenJet->pt() < 30.0 || std::abs(iterGenJet->eta()) > 2.5) continue;
        v_jets.push_back(iterGenJet->p4());
        v_jets_highpt.push_back((iterGenJet->pt() > 100.0) ? 1 : 0);
        v_jets_vhighpt.push_back((iterGenJet->pt() > 200.0) ? 1 : 0);
      }
    }
  }

  const int NJETS = 10;
  float jets_tp_sumpt[NJETS] = {0}; float jets_matchtrk_sumpt[NJETS] = {0}; float jets_trk_sumpt[NJETS] = {0};

  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------
  if (SaveAllTracks) {
    int this_l1track = 0;
    for (auto iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++) {
      edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> l1track_ptr(TTTrackHandle, this_l1track++);

      float tmp_trk_pt = iterL1Track->momentum().perp();
      float tmp_trk_eta = iterL1Track->momentum().eta();
      float tmp_trk_phi = iterL1Track->momentum().phi();
      float tmp_trk_z0 = iterL1Track->z0();  
      float tmp_trk_tanL = iterL1Track->tanL();
      int tmp_trk_charge = (int)TMath::Sign(1, iterL1Track->rInv());
      
      if (hphSetup->useNewKF() && fabs(tmp_trk_z0) > 30.) continue;

      int tmp_trk_hitpattern = (int)iterL1Track->hitPattern();
      hph::HitPatternHelper hph(hphSetup, tmp_trk_hitpattern, tmp_trk_tanL, tmp_trk_z0);
      std::vector<int> hp_binary = hph.binary();
      int tmp_trk_lhits_hitpattern = 0; int tmp_trk_dhits_hitpattern = 0;
      for (int i = 0; i < (int)hp_binary.size(); i++) {
        if (hp_binary[i]) {
          if (i < 6) tmp_trk_lhits_hitpattern += pow(10, i);
          else tmp_trk_dhits_hitpattern += pow(10, i - 6);
        }
      }

      float tmp_trk_d0 = (L1Tk_nPar == 5) ? (iterL1Track->POCA().x() * sin(tmp_trk_phi) - iterL1Track->POCA().y() * cos(tmp_trk_phi)) : 999.;
      std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>>> stubRefs = iterL1Track->getStubRefs();
      int tmp_trk_nstub = (int)stubRefs.size();
      
      int tmp_trk_dhits = 0; int tmp_trk_lhits = 0;
      for (int is = 0; is < tmp_trk_nstub; is++) {
        DetId detIdStub = theTrackerGeom->idToDet((stubRefs.at(is)->clusterRef(0))->getDetId())->geographicalId();
        int layer = static_cast<int>(tTopo->layer(detIdStub));
        if (detIdStub.subdetId() == StripSubdetector::TOB) tmp_trk_lhits += pow(10, layer - 1);
        else if (detIdStub.subdetId() == StripSubdetector::TID) tmp_trk_dhits += pow(10, layer - 1);
      }

      m_trk_pt->push_back(tmp_trk_pt); m_trk_eta->push_back(tmp_trk_eta); m_trk_phi->push_back(tmp_trk_phi);
      m_trk_z0->push_back(tmp_trk_z0); m_trk_d0->push_back(tmp_trk_d0); m_trk_chi2->push_back(iterL1Track->chi2());
      m_trk_chi2_dof->push_back(iterL1Track->chi2() / (2 * tmp_trk_nstub - L1Tk_nPar));
      m_trk_chi2rphi->push_back(iterL1Track->chi2XY()); m_trk_chi2rz->push_back(iterL1Track->chi2Z());
      m_trk_bendchi2->push_back(iterL1Track->stubPtConsistency()); m_trk_MVA1->push_back(iterL1Track->trkMVA1());
      m_trk_nstub->push_back(tmp_trk_nstub); m_trk_dhits->push_back(tmp_trk_dhits); m_trk_lhits->push_back(tmp_trk_lhits);
      m_trk_seed->push_back((int)iterL1Track->trackSeedType()); m_trk_hitpattern->push_back(tmp_trk_hitpattern);
      m_trk_lhits_hitpattern->push_back(tmp_trk_lhits_hitpattern); m_trk_dhits_hitpattern->push_back(tmp_trk_dhits_hitpattern);
      m_trk_nPSstub_hitpattern->push_back(hph.numPS()); m_trk_n2Sstub_hitpattern->push_back(hph.num2S());
      m_trk_nLostPSstub_hitpattern->push_back(hph.numMissingPS()); m_trk_nLost2Sstub_hitpattern->push_back(hph.numMissing2S());
      m_trk_nLoststub_V1_hitpattern->push_back(hph.numMissingInterior1()); m_trk_nLoststub_V2_hitpattern->push_back(hph.numMissingInterior2());
      m_trk_charge->push_back(tmp_trk_charge); m_trk_phiSector->push_back(iterL1Track->phiSector());
      m_trk_etaSector->push_back(hph.etaSector());
      m_trk_genuine->push_back(MCTruthTTTrackHandle->isGenuine(l1track_ptr) ? 1 : 0);
      m_trk_loose->push_back(MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr) ? 1 : 0);
      m_trk_unknown->push_back(MCTruthTTTrackHandle->isUnknown(l1track_ptr) ? 1 : 0);
      m_trk_combinatoric->push_back(MCTruthTTTrackHandle->isCombinatoric(l1track_ptr) ? 1 : 0);

      edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);
      int myFake = my_tp.isNull() ? 0 : ((my_tp->eventId().event() > 0) ? 2 : 1);
      m_trk_fake->push_back(myFake);

      float tmp_matchtp_pt = -999; float tmp_matchtp_eta = -999; float tmp_matchtp_phi = -999;
      if (!my_tp.isNull() && my_tp->eventId().event() == 0) {
        tmp_matchtp_pt = my_tp->pt(); tmp_matchtp_eta = my_tp->eta(); tmp_matchtp_phi = my_tp->phi();
      }
      m_trk_matchtp_pdgid->push_back(!my_tp.isNull() ? my_tp->pdgId() : -999);
      m_trk_matchtp_pt->push_back(tmp_matchtp_pt); m_trk_matchtp_eta->push_back(tmp_matchtp_eta); m_trk_matchtp_phi->push_back(tmp_matchtp_phi);
      m_trk_matchtp_z0->push_back(!my_tp.isNull() ? my_tp->vz() : -999); m_trk_matchtp_lxy->push_back(!my_tp.isNull() ? sqrt(my_tp->vx()*my_tp->vx() + my_tp->vy()*my_tp->vy()) : -999);
      m_trk_matchtp_d0->push_back(-999);

      if (TrackingInJets) {
        int InJet = 0; int InJetH = 0; int InJetVH = 0;
        for (int ij = 0; ij < (int)v_jets.size(); ij++) {
          float dR = deltaR(tmp_trk_eta, tmp_trk_phi, v_jets.at(ij).eta(), v_jets.at(ij).phi());
          if (dR < 0.4) {
            InJet = 1;
            if (v_jets_highpt.at(ij)) InJetH = 1;
            if (v_jets_vhighpt.at(ij)) InJetVH = 1;
            if (ij < NJETS) jets_trk_sumpt[ij] += tmp_trk_pt;
          }
        }
        m_trk_injet->push_back(InJet); m_trk_injet_highpt->push_back(InJetH); m_trk_injet_vhighpt->push_back(InJetVH);
      }

      const TTBV hitPattern((int)iterL1Track->hitPattern(), setup->numLayers());
      const double zT = iterL1Track->z0() + setup->chosenRofZ() * iterL1Track->tanL();
      const vector<int>& le = layerEncoding->layerEncoding(zT);
      vector<int> layers; layers.reserve(hitPattern.size());
      for (int layer : hitPattern.ids()) layers.push_back(le[layer]);
      m_trk_layers->push_back(layers);
    }
  }

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------
  int this_tp = 0;
  for (auto iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
    edm::Ptr<TrackingParticle> tp_ptr(TrackingParticleHandle, this_tp++);

    if (MyProcess != 1 && iterTP->eventId().event() > 0) continue;
    if (iterTP->pt() < TP_minPt || iterTP->charge() == 0. || std::abs(iterTP->eta()) > TP_maxEta) continue;

    int tmp_tp_pdgid = iterTP->pdgId();
    if (MyProcess == 13 && abs(tmp_tp_pdgid) != 13) continue;
    if (MyProcess == 11 && abs(tmp_tp_pdgid) != 11) continue;
    if ((MyProcess == 6 || MyProcess == 15 || MyProcess == 211) && abs(tmp_tp_pdgid) != 211) continue;

    float lxy = sqrt(iterTP->vx() * iterTP->vx() + iterTP->vy() * iterTP->vy());
    if (MyProcess == 6 && lxy > 1.0) continue;

    if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).empty()) continue;

    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>>> theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    int nStubTP = (int)theStubRefs.size();
    if (TP_minNStub > 0 && nStubTP < TP_minNStub) continue;

    int hasStubInLayer[11] = {0}; int nStubLayerTP = 0;
    for (auto& theStubRef : theStubRefs) {
      DetId detid(theStubRef->getDetId());
      int layer = (detid.subdetId() == StripSubdetector::TOB) ? static_cast<int>(tTopo->layer(detid)) - 1 : static_cast<int>(tTopo->layer(detid)) + 5;
      hasStubInLayer[layer] = MCTruthTTStubHandle->findTrackingParticlePtr(theStubRef).isNull() ? 1 : 2;
    }
    for (int isum : hasStubInLayer) if (isum >= 1) nStubLayerTP++;
    if (TP_minNStubLayer > 0 && nStubLayerTP < TP_minNStubLayer) continue;

    std::vector<edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>>> matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);
    int nMatch = 0; int i_track = -1; float i_chi2dof = 99999;

    for (int it = 0; it < (int)matchedTracks.size(); it++) {
      if (!MCTruthTTTrackHandle->isLooselyGenuine(matchedTracks.at(it))) continue;
      int tmp_trk_nstub = matchedTracks.at(it)->getStubRefs().size();
      if (tmp_trk_nstub < L1Tk_minNStub) continue;

      edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
      if (!my_tp.isNull() && std::abs(my_tp->pt() - iterTP->pt()) < 0.1 && tmp_tp_pdgid == my_tp->pdgId() && MCTruthTTTrackHandle->isGenuine(matchedTracks.at(it))) {
        nMatch++;
        float tmp_trk_chi2dof = matchedTracks.at(it)->chi2() / (2 * tmp_trk_nstub - L1Tk_nPar);
        if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
          i_track = it; i_chi2dof = tmp_trk_chi2dof;
        }
      }
    }

    m_tp_pt->push_back(iterTP->pt()); m_tp_eta->push_back(iterTP->eta()); m_tp_phi->push_back(iterTP->phi());
    m_tp_lxy->push_back(lxy); m_tp_z0->push_back(iterTP->vz()); m_tp_pdgid->push_back(tmp_tp_pdgid);
    m_tp_nmatch->push_back(nMatch); m_tp_nstub->push_back(nStubTP); m_tp_charge->push_back(iterTP->charge());

    if (nMatch > 0) {
      m_matchtrk_pt->push_back(matchedTracks.at(i_track)->momentum().perp());
      m_matchtrk_eta->push_back(matchedTracks.at(i_track)->momentum().eta());
      m_matchtrk_phi->push_back(matchedTracks.at(i_track)->momentum().phi());
      m_matchtrk_z0->push_back(matchedTracks.at(i_track)->z0());
      m_matchtrk_nstub->push_back((int)matchedTracks.at(i_track)->getStubRefs().size());
    } else {
      m_matchtrk_pt->push_back(-999); m_matchtrk_eta->push_back(-999); m_matchtrk_phi->push_back(-999);
      m_matchtrk_z0->push_back(-999); m_matchtrk_nstub->push_back(-999);
    }
  }

  if (TrackingInJets) {
    for (int ij = 0; ij < (int)v_jets.size() && ij < NJETS; ij++) {
      m_jet_eta->push_back(v_jets.at(ij).eta()); m_jet_phi->push_back(v_jets.at(ij).phi()); m_jet_pt->push_back(v_jets.at(ij).pt());
      m_jet_tp_sumpt->push_back(jets_tp_sumpt[ij]); m_jet_trk_sumpt->push_back(jets_trk_sumpt[ij]);
    }
  }

  eventTree->Fill();
}

DEFINE_FWK_MODULE(L1TrackHitNtupleMaker);

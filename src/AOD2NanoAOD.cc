// -*- C++ -*-
//
// Package:    AOD2NanoAOD
// Class:      AOD2NanoAOD

#include <memory>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class AOD2NanoAOD : public edm::EDAnalyzer {
public:
  explicit AOD2NanoAOD(const edm::ParameterSet &);
  ~AOD2NanoAOD();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  bool providesGoodLumisection(const edm::Event &iEvent);

  TTree *tree;

  // Event information
  Int_t value_run;
  UInt_t value_lumi_block;
  ULong64_t value_event;

  // Vertices
  int value_ve_n;
  float value_ve_x;
  float value_ve_y;
  float value_ve_z;

  // Muons
  const static int max_mu = 100;
  UInt_t value_mu_n;
  float value_mu_pt[max_mu];
  float value_mu_eta[max_mu];
  float value_mu_phi[max_mu];
  float value_mu_mass[max_mu];
  int value_mu_charge[max_mu];

  // Electrons
  const static int max_el = 100;
  UInt_t value_el_n;
  float value_el_pt[max_el];
  float value_el_eta[max_el];
  float value_el_phi[max_el];
  float value_el_mass[max_el];
  int value_el_charge[max_el];

  // Taus
  const static int max_tau = 300;
  UInt_t value_tau_n;
  float value_tau_pt[max_tau];
  float value_tau_eta[max_tau];
  float value_tau_phi[max_tau];
  float value_tau_mass[max_tau];
  int value_tau_charge[max_tau];

  // Photons
  const static int max_ph = 100;
  UInt_t value_ph_n;
  float value_ph_pt[max_ph];
  float value_ph_eta[max_ph];
  float value_ph_phi[max_ph];
  float value_ph_mass[max_ph];
  int value_ph_charge[max_ph];

  // MET
  float value_met_pt;
  float value_met_phi;
  float value_met_sumet;
  float value_met_significance;
  float value_met_covxx;
  float value_met_covxy;
  float value_met_covyy;

  // Jets
  const static int max_jet = 300;
  UInt_t value_jet_n;
  float value_jet_pt[max_jet];
  float value_jet_eta[max_jet];
  float value_jet_phi[max_jet];
  float value_jet_mass[max_jet];
};

AOD2NanoAOD::AOD2NanoAOD(const edm::ParameterSet &iConfig) {
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("Events", "Events");

  // Event information
  tree->Branch("run", &value_run);
  tree->Branch("luminosityBlock", &value_lumi_block);
  tree->Branch("event", &value_event);

  // Vertices
  tree->Branch("PV_npvs", &value_ve_n, "PV_npvs/I");
  tree->Branch("PV_x", &value_ve_x, "PV_x/F");
  tree->Branch("PV_y", &value_ve_y, "PV_y/F");
  tree->Branch("PV_z", &value_ve_z, "PV_z/F");

  // Muons
  tree->Branch("nMuon", &value_mu_n, "nMuon/i");
  tree->Branch("Muon_pt", value_mu_pt, "Muon_pt[nMuon]/F");
  tree->Branch("Muon_eta", value_mu_eta, "Muon_eta[nMuon]/F");
  tree->Branch("Muon_phi", value_mu_phi, "Muon_phi[nMuon]/F");
  tree->Branch("Muon_mass", value_mu_mass, "Muon_mass[nMuon]/F");
  tree->Branch("Muon_charge", value_mu_charge, "Muon_charge[nMuon]/I");

  // Electrons
  tree->Branch("nElectron", &value_el_n, "nElectron/i");
  tree->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_mass", value_el_mass, "Electron_mass[nElectron]/F");
  tree->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");

  // Taus
  tree->Branch("nTau", &value_tau_n, "nTau/i");
  tree->Branch("Tau_pt", value_tau_pt, "Tau_pt[nTau]/F");
  tree->Branch("Tau_eta", value_tau_eta, "Tau_eta[nTau]/F");
  tree->Branch("Tau_phi", value_tau_phi, "Tau_phi[nTau]/F");
  tree->Branch("Tau_mass", value_tau_mass, "Tau_mass[nTau]/F");
  tree->Branch("Tau_charge", value_tau_charge, "Tau_charge[nTau]/I");

  // Photons
  tree->Branch("nPhoton", &value_ph_n, "nPhoton/i");
  tree->Branch("Photon_pt", value_ph_pt, "Photon_pt[nPhoton]/F");
  tree->Branch("Photon_eta", value_ph_eta, "Photon_eta[nPhoton]/F");
  tree->Branch("Photon_phi", value_ph_phi, "Photon_phi[nPhoton]/F");
  tree->Branch("Photon_mass", value_ph_mass, "Photon_mass[nPhoton]/F");
  tree->Branch("Photon_charge", value_ph_charge, "Photon_charge[nPhoton]/I");

  // MET
  tree->Branch("MET_pt", &value_met_pt, "MET_pt/F");
  tree->Branch("MET_phi", &value_met_phi, "MET_phi/F");
  tree->Branch("MET_sumet", &value_met_sumet, "MET_sumet/F");
  tree->Branch("MET_significance", &value_met_significance, "MET_significance/F");
  tree->Branch("MET_CovXX", &value_met_covxx, "MET_CovXX/F");
  tree->Branch("MET_CovXY", &value_met_covxy, "MET_CovXY/F");
  tree->Branch("MET_CovYY", &value_met_covyy, "MET_CovYY/F");

  // Jets
  tree->Branch("nJet", &value_jet_n, "nJet/i");
  tree->Branch("Jet_pt", value_jet_pt, "Jet_pt[nJet]/F");
  tree->Branch("Jet_eta", value_jet_eta, "Jet_eta[nJet]/F");
  tree->Branch("Jet_phi", value_jet_phi, "Jet_phi[nJet]/F");
  tree->Branch("Jet_mass", value_jet_mass, "Jet_mass[nJet]/F");
}

AOD2NanoAOD::~AOD2NanoAOD() {}

void AOD2NanoAOD::analyze(const edm::Event &iEvent,
                          const edm::EventSetup &iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  /*
  LogInfo("DEBUG")
  << "Starting to analyze \n"
  << "Event number: " << (iEvent.id()).event()
  << ", Run number: " << iEvent.run()
  << ", Lumisection: " << iEvent.luminosityBlock();
  */


  // Event information
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#GetData
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#get_ByLabel
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable
  value_run = iEvent.run();
  value_lumi_block = iEvent.luminosityBlock();
  value_event = iEvent.id().event();

  // Vertex
  Handle<VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  value_ve_n = vertices->size();
  value_ve_x = vertices->begin()->x();
  value_ve_y = vertices->begin()->y();
  value_ve_z = vertices->begin()->z();

  // Muons
  Handle<MuonCollection> muons;
  iEvent.getByLabel(InputTag("muons"), muons);

  value_mu_n = muons->size();
  for (auto it = muons->begin(); it != muons->end(); it++) {
    const auto idx = it - muons->begin();
    value_mu_pt[idx] = it->pt();
    value_mu_eta[idx] = it->eta();
    value_mu_phi[idx] = it->phi();
    value_mu_charge[idx] = it->charge();
    value_mu_mass[idx] = it->mass();
  }

  // Electrons
  Handle<GsfElectronCollection> electrons;
  iEvent.getByLabel(InputTag("gsfElectrons"), electrons);

  value_el_n = electrons->size();
  for (auto it = electrons->begin(); it != electrons->end(); it++) {
    const auto idx = it - electrons->begin();
    value_el_pt[idx] = it->pt();
    value_el_eta[idx] = it->eta();
    value_el_phi[idx] = it->phi();
    value_el_charge[idx] = it->charge();
    value_el_mass[idx] = it->mass();
  }

  // Taus
  Handle<PFTauCollection> taus;
  iEvent.getByLabel(InputTag("hpsPFTauProducer"), taus);

  value_tau_n = taus->size();
  for (auto it = taus->begin(); it != taus->end(); it++) {
    const auto idx = it - taus->begin();
    value_tau_pt[idx] = it->pt();
    value_tau_eta[idx] = it->eta();
    value_tau_phi[idx] = it->phi();
    value_tau_charge[idx] = it->charge();
    value_tau_mass[idx] = it->mass();
  }

  // Photons
  Handle<PhotonCollection> photons;
  iEvent.getByLabel(InputTag("photons"), photons);

  value_ph_n = photons->size();
  for (auto it = photons->begin(); it != photons->end(); it++) {
    const auto idx = it - photons->begin();
    value_ph_pt[idx] = it->pt();
    value_ph_eta[idx] = it->eta();
    value_ph_phi[idx] = it->phi();
    value_ph_charge[idx] = it->charge();
    value_ph_mass[idx] = it->mass();
  }

  // MET
  Handle<PFMETCollection> met;
  iEvent.getByLabel(InputTag("pfMet"), met);
  value_met_pt = met->begin()->pt();
  value_met_phi = met->begin()->phi();
  value_met_sumet = met->begin()->sumEt();
  value_met_significance = met->begin()->significance();
  auto cov = met->begin()->getSignificanceMatrix();
  value_met_covxx = cov[0][0];
  value_met_covxy = cov[0][1];
  value_met_covyy = cov[1][1];

  // Jets
  Handle<PFJetCollection> jets;
  iEvent.getByLabel(InputTag("ak5PFJets"), jets);

  value_jet_n = jets->size();
  for (auto it = jets->begin(); it != jets->end(); it++) {
    const auto idx = it - jets->begin();
    value_jet_pt[idx] = it->pt();
    value_jet_eta[idx] = it->eta();
    value_jet_phi[idx] = it->phi();
    value_jet_mass[idx] = it->mass();
  }

  // Fill event
  tree->Fill();
}

void AOD2NanoAOD::beginJob() {}

void AOD2NanoAOD::endJob() {}

DEFINE_FWK_MODULE(AOD2NanoAOD);

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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/ParameterSet/interface/Registry.h"

const static std::vector<std::string> interestingTriggers = {
        "HLT_Activity_Ecal_SC7_v12",
        "HLT_L1SingleJet16_v6",
        "HLT_L1SingleJet36_v6",
        "HLT_Jet20_NoL1FastJet_v2",
        "HLT_PFJet40_v5",
        "HLT_Jet50_NoL1FastJet_v2",
        "HLT_PFJet80_v5",
        "HLT_PFJet140_v5",
        "HLT_PFJet200_v5",
        "HLT_PFJet260_v5",
        "HLT_PFJet320_v5",
        "HLT_Jet370_NoJetID_v13",
        "HLT_PFJet400_v5",
        "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5",
        "HLT_SingleJetC5_v2",
        "HLT_SingleForJet25_v2",
        "HLT_SingleForJet15_v2",
        "HLT_DiPFJetAve40_v6",
        "HLT_DiPFJetAve80_v6",
        "HLT_DiPFJetAve140_v6",
        "HLT_DiPFJetAve200_v6",
        "HLT_DiPFJetAve260_v6",
        "HLT_DiPFJetAve320_v6",
        "HLT_DiPFJetAve400_v6",
        "HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5",
        "HLT_DoubleJet20_ForwardBackward_v2",
        "HLT_DiJet80_DiJet60_DiJet20_v2",
        "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v5",
        "HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v5",
        "HLT_DiJet40Eta2p6_BTagIP3DFastPV_v4",
        "HLT_DiJet80Eta2p6_BTagIP3DFastPVLoose_v4",
        "HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v4",
        "HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV_v4",
        "HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose_v4",
        "HLT_QuadJet50_v2",
        "HLT_QuadJet50_Jet20_v1",
        "HLT_QuadJet60_DiJet20_v2",
        "HLT_QuadJet70_v2",
        "HLT_QuadJet80_v2",
        "HLT_QuadJet90_v2",
        "HLT_QuadJet75_55_35_20_BTagIP_VBF_v3",
        "HLT_QuadJet75_55_38_20_BTagIP_VBF_v3",
        "HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v1",
        "HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v1",
        "HLT_SixJet35_v2",
        "HLT_SixJet45_v2",
        "HLT_SixJet50_v2",
        "HLT_EightJet30_eta3p0_v1",
        "HLT_EightJet35_eta3p0_v1",
        "HLT_ExclDiJet35_HFOR_v2",
        "HLT_ExclDiJet35_HFAND_v2",
        "HLT_ExclDiJet80_HFAND_v2",
        "HLT_JetE30_NoBPTX_v12",
        "HLT_JetE30_NoBPTX3BX_v1",
        "HLT_JetE50_NoBPTX3BX_v1",
        "HLT_JetE70_NoBPTX3BX_v1",
        "HLT_HT200_AlphaT0p57_v5",
        "HLT_HT200_v3",
        "HLT_HT250_AlphaT0p55_v4",
        "HLT_HT250_AlphaT0p57_v4",
        "HLT_HT250_v3",
        "HLT_HT300_AlphaT0p53_v4",
        "HLT_HT300_AlphaT0p54_v10",
        "HLT_HT300_v3",
        "HLT_HT300_DoubleDisplacedPFJet60_v5",
        "HLT_HT300_DoubleDisplacedPFJet60_ChgFraction10_v5",
        "HLT_HT300_SingleDisplacedPFJet60_v5",
        "HLT_HT300_SingleDisplacedPFJet60_ChgFraction10_v5",
        "HLT_HT350_v3",
        "HLT_HT350_AlphaT0p52_v4",
        "HLT_HT350_AlphaT0p53_v15",
        "HLT_HT400_v3",
        "HLT_HT400_AlphaT0p51_v15",
        "HLT_HT400_AlphaT0p52_v10",
        "HLT_HT450_AlphaT0p51_v10",
        "HLT_HT450_v3",
        "HLT_HT500_v3",
        "HLT_HT550_v3",
        "HLT_HT650_v3",
        "HLT_HT650_Track50_dEdx3p6_v6",
        "HLT_HT650_Track60_dEdx3p7_v6",
        "HLT_HT750_v3",
        "HLT_PFHT350_v6",
        "HLT_PFHT650_v8",
        "HLT_PFHT650_DiCentralPFJet80_CenPFJet40_v6",
        "HLT_PFHT700_v6",
        "HLT_PFHT750_v6",
        "HLT_PFMET150_v4",
        "HLT_PFMET180_v4",
        "HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v1",
        "HLT_DiCentralPFJet30_PFMET80_v2",
        "HLT_DiCentralPFJet50_PFMET80_v6",
        "HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v1",
        "HLT_DiPFJet80_DiPFJet30_BTagCSVd07d05d03_v1",
        "HLT_DiPFJet80_DiPFJet30_BTagCSVd07d05d05_v1",
        "HLT_DiPFJet80_DiPFJet30_BTagCSVd07d05d03_PFDiJetPt120_v1",
        "HLT_DiPFJet80_DiPFJet30_BTagCSVd07d05_v1",
        "HLT_MET80_v3",
        "HLT_MET80_Track50_dEdx3p6_v4",
        "HLT_MET80_Track60_dEdx3p7_v4",
        "HLT_MET120_v10",
        "HLT_MET120_HBHENoiseCleaned_v3",
        "HLT_MET200_v10",
        "HLT_MET200_HBHENoiseCleaned_v3",
        "HLT_MET300_v2",
        "HLT_MET300_HBHENoiseCleaned_v3",
        "HLT_MET400_v5",
        "HLT_MET400_HBHENoiseCleaned_v3",
        "HLT_L1SingleMuOpen_v6",
        "HLT_L1SingleMu12_v1",
        "HLT_L2Mu70_eta2p1_PFMET55_v1",
        "HLT_L2Mu70_eta2p1_PFMET60_v1",
        "HLT_L2Mu20_eta2p1_NoVertex_v1",
        "HLT_L2Mu10_NoVertex_NoBPTX3BX_v1",
        "HLT_L2Mu20_NoVertex_NoBPTX3BX_v1",
        "HLT_L2Mu30_NoVertex_NoBPTX3BX_v1",
        "HLT_L2TripleMu10_0_0_NoVertex_PFJet40Neutral_v4",
        "HLT_DoubleDisplacedMu4_DiPFJet40Neutral_v4",
        "HLT_Mu5_v18",
        "HLT_Mu8_v16",
        "HLT_Mu12_v16",
        "HLT_Mu17_v3",
        "HLT_Mu12_eta2p1_L1Mu10erJetC12WdEtaPhi1DiJetsC_v3",
        "HLT_Mu15_eta2p1_v3",
        "HLT_Mu24_v14",
        "HLT_Mu24_eta2p1_v3",
        "HLT_Mu30_v14",
        "HLT_Mu30_eta2p1_v3",
        "HLT_Mu40_v12",
        "HLT_Mu40_eta2p1_v9",
        "HLT_Mu50_eta2p1_v6",
        "HLT_RelIso1p0Mu5_v4",
        "HLT_RelIso1p0Mu20_v1",
        "HLT_IsoMu15_eta2p1_L1ETM20_v5",
        "HLT_IsoMu20_eta2p1_v5",
        "HLT_IsoMu24_v15",
        "HLT_IsoMu24_eta2p1_v13",
        "HLT_IsoMu30_v9",
        "HLT_IsoMu30_eta2p1_v13",
        "HLT_IsoMu34_eta2p1_v11",
        "HLT_IsoMu40_eta2p1_v8",
        "HLT_Mu40_eta2p1_Track50_dEdx3p6_v3",
        "HLT_Mu40_eta2p1_Track60_dEdx3p7_v3",
        "HLT_L2DoubleMu23_NoVertex_v10",
        "HLT_L2DoubleMu23_NoVertex_2Cha_Angle2p5_v2",
        "HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_v2",
        "HLT_DoubleMu11_Acoplanarity03_v3",
        "HLT_DoubleMu4_Jpsi_Displaced_v9",
        "HLT_DoubleMu4_JpsiTk_Displaced_v3",
        "HLT_DoubleMu3_4_Dimuon5_Bs_Central_v2",
        "HLT_DoubleMu3p5_4_Dimuon5_Bs_Central_v2",
        "HLT_DoubleMu4_Dimuon7_Bs_Forward_v2",
        "HLT_DoubleMu3p5_LowMass_Displaced_v3",
        "HLT_DoubleMu3p5_LowMassNonResonant_Displaced_v3",
        "HLT_Dimuon0_Jpsi_v14",
        "HLT_Dimuon0_Jpsi_NoVertexing_v11",
        "HLT_Dimuon0_Upsilon_v14",
        "HLT_Dimuon0_PsiPrime_v3",
        "HLT_Dimuon5_Upsilon_v3",
        "HLT_Dimuon5_PsiPrime_v3",
        "HLT_Dimuon7_Upsilon_v4",
        "HLT_Dimuon8_Jpsi_v4",
        "HLT_Dimuon8_Upsilon_v3",
        "HLT_Dimuon9_PsiPrime_v9",
        "HLT_Dimuon10_Jpsi_v3",
        "HLT_Dimuon11_Upsilon_v3",
        "HLT_Dimuon0_Jpsi_Muon_v15",
        "HLT_Dimuon0_Upsilon_Muon_v15",
        "HLT_Dimuon3p5_SameSign_v3",
        "HLT_DoubleMu4_Acoplanarity03_v3",
        "HLT_Tau2Mu_ItTrack_v3",
        "HLT_Mu13_Mu8_v17",
        "HLT_Mu17_Mu8_v17",
        "HLT_Mu17_TkMu8_v10",
        "HLT_Mu22_TkMu8_v6",
        "HLT_Mu22_TkMu22_v6",
        "HLT_TripleMu5_v17",
        "HLT_DoubleMu5_IsoMu5_v18",
        "HLT_Mu5_L2Mu3_Jpsi_v4",
        "HLT_Mu5_Track2_Jpsi_v18",
        "HLT_Mu5_Track3p5_Jpsi_v4",
        "HLT_Mu7_Track7_Jpsi_v18",
        "HLT_Photon20_CaloIdVL_v3",
        "HLT_Photon20_CaloIdVL_IsoL_v15",
        "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_v4",
        "HLT_Photon26_Photon18_v11",
        "HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v3",
        "HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v5",
        "HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v5",
        "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v5",
        "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v5",
        "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v1",
        "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v4",
        "HLT_Photon30_CaloIdVL_v13",
        "HLT_Photon30_CaloIdVL_IsoL_v18",
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v4",
        "HLT_Photon36_Photon22_v5",
        "HLT_Photon36_R9Id85_Photon22_R9Id85_v3",
        "HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v5",
        "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v5",
        "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v5",
        "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v5",
        "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v4",
        "HLT_Photon50_CaloIdVL_v9",
        "HLT_Photon50_CaloIdVL_IsoL_v16",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v4",
        "HLT_Photon60_CaloIdL_MHT70_v8",
        "HLT_Photon60_CaloIdL_HT300_v1",
        "HLT_Photon70_CaloIdXL_PFHT400_v5",
        "HLT_Photon70_CaloIdXL_PFHT500_v5",
        "HLT_Photon70_CaloIdXL_PFMET100_v4",
        "HLT_Photon75_CaloIdVL_v12",
        "HLT_Photon75_CaloIdVL_IsoL_v17",
        "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v4",
        "HLT_Photon90_CaloIdVL_v9",
        "HLT_Photon90_CaloIdVL_IsoL_v14",
        "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v4",
        "HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25_v1",
        "HLT_DisplacedPhoton65EBOnly_CaloIdVL_IsoL_PFMET30_v1",
        "HLT_Photon135_v6",
        "HLT_Photon150_v3",
        "HLT_Photon160_v3",
        "HLT_Photon300_NoHE_v4",
        "HLT_DoublePhoton48_HEVT_v7",
        "HLT_DoublePhoton53_HEVT_v1",
        "HLT_DoublePhoton70_v5",
        "HLT_DoublePhoton80_v6",
        "HLT_DoublePhoton5_IsoVL_CEP_v15",
        "HLT_L1SingleEG5_v5",
        "HLT_L1SingleEG12_v5",
        "HLT_L1DoubleEG3_FwdVeto_v1",
        "HLT_L1ETM30_v1",
        "HLT_L1ETM40_v1",
        "HLT_L1ETM70_v1",
        "HLT_L1ETM100_v1",
        "HLT_Ele8_CaloIdT_TrkIdVL_v4",
        "HLT_Ele8_CaloIdT_TrkIdVL_EG7_v1",
        "HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v4",
        "HLT_Ele8_CaloIdL_CaloIsoVL_v16",
        "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14",
        "HLT_Ele17_CaloIdL_CaloIsoVL_v16",
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",
        "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v5",
        "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v5",
        "HLT_Ele22_CaloIdL_CaloIsoVL_v5",
        "HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",
        "HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele15_CaloIdT_CaloIsoVL_trackless_v6",
        "HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT15_v6",
        "HLT_Ele23_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT30_v6",
        "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet30_v1",
        "HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet50_40_30_v1",
        "HLT_Ele27_WP80_v10",
        "HLT_Ele27_WP80_PFMET_MT50_v4",
        "HLT_Ele30_CaloIdVT_TrkIdT_v5",
        "HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",
        "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5",
        "HLT_Ele80_CaloIdVT_GsfTrkIdT_v1",
        "HLT_Ele90_CaloIdVT_GsfTrkIdT_v1",
        "HLT_DoubleEle8_CaloIdT_TrkIdVL_v11",
        "HLT_DoubleEle33_CaloIdL_v13",
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6",
        "HLT_DoubleEle33_CaloIdT_v9",
        "HLT_Ele15_Ele8_Ele5_CaloIdL_TrkIdVL_v5",
        "HLT_LooseIsoPFTau35_Trk20_Prong1_v6",
        "HLT_LooseIsoPFTau35_Trk20_Prong1_MET70_v6",
        "HLT_LooseIsoPFTau35_Trk20_Prong1_MET75_v6",
        "HLT_IsoMu15_eta2p1_LooseIsoPFTau35_Trk20_Prong1_L1ETM20_v6",
        "HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v2",
        "HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_v2",
        "HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v6",
        "HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1_v6",
        "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v7",
        "HLT_BTagMu_DiJet20_Mu5_v3",
        "HLT_BTagMu_DiJet40_Mu5_v3",
        "HLT_BTagMu_DiJet70_Mu5_v3",
        "HLT_BTagMu_DiJet110_Mu5_v3",
        "HLT_BTagMu_Jet300_Mu5_v3",
        "HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v5",
        "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",
        "HLT_Mu8_DiJet30_v4",
        "HLT_Mu8_TriJet30_v4",
        "HLT_Mu8_QuadJet30_v4",
        "HLT_IsoMu12_DoubleCentralJet65_v1",
        "HLT_Mu15_eta2p1_L1ETM20_v3",
        "HLT_Mu24_PFJet30_PFJet25_Deta3_CentralPFJet25_v1",
        "HLT_Mu24_CentralPFJet30_CentralPFJet25_v1",
        "HLT_IsoMu24_PFJet30_PFJet25_Deta3_CentralPFJet25_v1",
        "HLT_IsoMu24_CentralPFJet30_CentralPFJet25_v1",
        "HLT_IsoMu24_CentralPFJet30_CentralPFJet25_PFMET20_v1",
        "HLT_Ele32_WP80_PFJet30_PFJet25_Deta3_v1",
        "HLT_Ele32_WP80_PFJet30_PFJet25_Deta3_CentralPFJet30_v1",
        "HLT_Ele32_WP80_CentralPFJet35_CentralPFJet25_v1",
        "HLT_Ele32_WP80_CentralPFJet35_CentralPFJet25_PFMET20_v1",
        "HLT_IsoMu17_eta2p1_CentralPFNoPUJet30_BTagIPIter_v1",
        "HLT_IsoMu17_eta2p1_CentralPFNoPUJet30_v1",
        "HLT_IsoMu17_eta2p1_DiCentralPFNoPUJet30_v1",
        "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1",
        "HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet50_40_30_v1",
        "HLT_Mu17_eta2p1_CentralPFNoPUJet30_BTagIPIter_v1",
        "HLT_Mu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1",
        "HLT_Mu17_eta2p1_TriCentralPFNoPUJet50_40_30_v1",
        "HLT_IsoMu20_WCandPt80_v1",
        "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",
        "HLT_Mu12_eta2p1_DiCentral_40_20_DiBTagIP3D1stTrack_v3",
        "HLT_Mu12_eta2p1_DiCentral_40_20_BTagIP3D1stTrack_v3",
        "HLT_Mu12_eta2p1_DiCentral_40_20_v3",
        "HLT_Mu12_eta2p1_DiCentral_20_v3",
        "HLT_Mu15_eta2p1_TriCentral_40_20_20_DiBTagIP3D1stTrack_v3",
        "HLT_Mu15_eta2p1_TriCentral_40_20_20_BTagIP3D1stTrack_v3",
        "HLT_Mu15_eta2p1_TriCentral_40_20_20_v3",
        "HLT_Mu30_Ele30_CaloIdL_v6",
        "HLT_IsoMu17_eta2p1_DiCentralPFJet30_PFHT350_PFMHT40_v6",
        "HLT_IsoMu20_eta2p1_CentralPFJet80_v6",
        "HLT_DoubleRelIso1p0Mu5_Mass8_PFHT175_v6",
        "HLT_DoubleRelIso1p0Mu5_Mass8_PFHT225_v6",
        "HLT_DoubleMu8_Mass8_PFHT225_v6",
        "HLT_DoubleMu8_Mass8_PFHT175_v6",
        "HLT_RelIso1p0Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT175_v6",
        "HLT_RelIso1p0Mu5_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT225_v6",
        "HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT175_v6",
        "HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass8_PFHT225_v6",
        "HLT_PFHT350_Mu15_PFMET45_v6",
        "HLT_PFHT350_Mu15_PFMET50_v6",
        "HLT_PFHT400_Mu5_PFMET45_v6",
        "HLT_PFHT400_Mu5_PFMET50_v6",
        "HLT_Mu40_PFHT350_v6",
        "HLT_Mu60_PFHT350_v6",
        "HLT_Mu40_HT200_v1",
        "HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v14",
        "HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v3",
        "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5",
        "HLT_Ele8_CaloIdT_TrkIdT_DiJet30_v15",
        "HLT_Ele8_CaloIdT_TrkIdT_TriJet30_v15",
        "HLT_Ele8_CaloIdT_TrkIdT_QuadJet30_v15",
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5",
        "HLT_Ele25_CaloIdVT_TrkIdT_TriCentralPFNoPUJet30_30_20_v1",
        "HLT_Ele25_CaloIdVT_TrkIdT_TriCentralPFNoPUJet50_40_30_v5",
        "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_v5",
        "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralPFNoPUJet30_v5",
        "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet30_30_20_v1",
        "HLT_Ele25_CaloIdVL_CaloIsoT_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet30_30_20_v1",
        "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_BTagIPIter_v6",
        "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet50_40_30_v5",
        "HLT_Ele27_WP80_CentralPFJet80_v6",
        "HLT_Ele27_WP80_WCandPt80_v6",
        "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v5",
        "HLT_Ele30_CaloIdVT_TrkIdT_PFJet150_PFJet25_v5",
        "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v5",
        "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet150_PFNoPUJet25_v5",
        "HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT175_v6",
        "HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass8_PFHT225_v6",
        "HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_CaloIdT_TrkIdVL_v11",
        "HLT_TripleEle10_CaloIdL_TrkIdVL_v17",
        "HLT_RsqMR40_Rsq0p04_v3",
        "HLT_RsqMR45_Rsq0p09_v2",
        "HLT_RsqMR55_Rsq0p09_MR150_v3",
        "HLT_RsqMR60_Rsq0p09_MR150_v3",
        "HLT_RsqMR65_Rsq0p09_MR150_v2",
        "HLT_IsoMu12_RsqMR30_Rsq0p04_MR200_v1",
        "HLT_IsoMu12_RsqMR40_Rsq0p04_MR200_v1",
        "HLT_Ele12_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_RsqMR30_Rsq0p04_MR200_v1",
        "HLT_Ele12_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_RsqMR40_Rsq0p04_MR200_v1",
        "HLT_Ele12_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_DoubleCentralJet65_v1",
        "HLT_Photon40_CaloIdL_RsqMR35_Rsq0p09_MR150_v3",
        "HLT_Photon40_CaloIdL_RsqMR40_Rsq0p09_MR150_v3",
        "HLT_Photon40_CaloIdL_RsqMR45_Rsq0p09_MR150_v3",
        "HLT_Photon40_CaloIdL_RsqMR50_Rsq0p09_MR150_v3",
        "HLT_DoublePhoton40_CaloIdL_Rsq0p035_v3",
        "HLT_DoublePhoton40_CaloIdL_Rsq0p06_v3",
        "HLT_Mu22_Photon22_CaloIdL_v5",
        "HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v5",
        "HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v5",
        "HLT_DoubleMu14_Mass8_PFMET40_v5",
        "HLT_DoubleMu14_Mass8_PFMET50_v5",
        "HLT_DoubleEle14_CaloIdT_TrkIdVL_Mass8_PFMET40_v5",
        "HLT_DoubleEle14_CaloIdT_TrkIdVL_Mass8_PFMET50_v5",
        "HLT_Mu14_Ele14_CaloIdT_TrkIdVL_Mass8_PFMET40_v5",
        "HLT_Mu14_Ele14_CaloIdT_TrkIdVL_Mass8_PFMET50_v5",
        "HLT_PFHT350_PFMET100_v6",
        "HLT_PFHT400_PFMET100_v6",
        "HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v5",
        "HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v5",
        "HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v5",
        "HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v5",
        "HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v5",
        "HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v5",
        "HLT_Ele5_SC5_Jpsi_Mass2to15_v3",
        "HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF_v1",
        "HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v1",
        "HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF_v1",
        "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v2",
        "HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_v2",
        "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2",
        "HLT_Mu17_eta2p1_LooseIsoPFTau20_v2",
        "HLT_PixelTracks_Multiplicity70_v2",
        "HLT_PixelTracks_Multiplicity80_v11",
        "HLT_PixelTracks_Multiplicity90_v2",
        "HLT_BeamGas_HF_Beam1_v4",
        "HLT_BeamGas_HF_Beam2_v4",
        "HLT_BeamHalo_v12",
        "HLT_IsoTrackHE_v14",
        "HLT_IsoTrackHB_v13",
        "HLT_HcalPhiSym_v10",
        "HLT_HcalNZS_v9",
        "HLT_GlobalRunHPDNoise_v7",
        "HLT_L1Tech_HBHEHO_totalOR_v5",
        "HLT_L1Tech_HCAL_HF_single_channel_v3",
        "HLT_ZeroBias_v6",
        "HLT_ZeroBiasPixel_DoubleTrack_v1",
        "HLT_Physics_v4",
        "HLT_DTCalibration_v2",
        "HLT_EcalCalibration_v3",
        "HLT_HcalCalibration_v3",
        "HLT_TrackerCalibration_v3",
        "HLT_Random_v2",
        "HLT_L1SingleMuOpen_AntiBPTX_v6",
        "HLT_L1TrackerCosmics_v6",
        "HLT_LogMonitor_v3",
        "HLT_DTErrors_v3",
        "HLT_L1DoubleJet36Central_v6",
};

class AOD2NanoAOD : public edm::EDAnalyzer {
public:
  explicit AOD2NanoAOD(const edm::ParameterSet &);
  ~AOD2NanoAOD();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  bool providesGoodLumisection(const edm::Event &iEvent);
  bool isData;

  TTree *tree;

  // Event information
  Int_t value_run;
  UInt_t value_lumi_block;
  ULong64_t value_event;

  // Trigger
  const static int max_trig = 1000;
  bool value_trig[max_trig];

  // Vertices
  int value_ve_n;
  float value_ve_x;
  float value_ve_y;
  float value_ve_z;

  // Muons
  const static int max_mu = 1000;
  UInt_t value_mu_n;
  float value_mu_pt[max_mu];
  float value_mu_eta[max_mu];
  float value_mu_phi[max_mu];
  float value_mu_mass[max_mu];
  int value_mu_charge[max_mu];
  float value_mu_pfreliso03all[max_mu];
  float value_mu_pfreliso04all[max_mu];
  bool value_mu_tightid[max_mu];
  bool value_mu_softid[max_mu];
  float value_mu_dxy[max_mu];
  float value_mu_dxyErr[max_mu];
  float value_mu_dz[max_mu];
  float value_mu_dzErr[max_mu];
  int value_mu_genpartidx[max_mu];
  int value_mu_jetidx[max_mu];

  // Electrons
  const static int max_el = 1000;
  UInt_t value_el_n;
  float value_el_pt[max_el];
  float value_el_eta[max_el];
  float value_el_phi[max_el];
  float value_el_mass[max_el];
  int value_el_charge[max_el];
  float value_el_pfreliso03all[max_el];
  float value_el_dxy[max_el];
  float value_el_dxyErr[max_el];
  float value_el_dz[max_el];
  float value_el_dzErr[max_el];
  bool value_el_cutbasedid[max_el];
  bool value_el_pfid[max_el];
  int value_el_genpartidx[max_el];
  int value_el_jetidx[max_el];

  // Taus
  const static int max_tau = 1000;
  UInt_t value_tau_n;
  float value_tau_pt[max_tau];
  float value_tau_eta[max_tau];
  float value_tau_phi[max_tau];
  float value_tau_mass[max_tau];
  int value_tau_charge[max_tau];
  int value_tau_decaymode[max_tau];
  float value_tau_chargediso[max_tau];
  float value_tau_neutraliso[max_tau];
  float value_tau_reliso_all[max_tau];
  int value_tau_genpartidx[max_tau];
  int value_tau_jetidx[max_tau];

  // Photons
  const static int max_ph = 1000;
  UInt_t value_ph_n;
  float value_ph_pt[max_ph];
  float value_ph_eta[max_ph];
  float value_ph_phi[max_ph];
  float value_ph_mass[max_ph];
  int value_ph_charge[max_ph];
  float value_ph_pfreliso03all[max_ph];
  int value_ph_genpartidx[max_ph];
  int value_ph_jetidx[max_ph];

  // MET
  float value_met_pt;
  float value_met_phi;
  float value_met_sumet;
  float value_met_significance;
  float value_met_covxx;
  float value_met_covxy;
  float value_met_covyy;

  // Jets
  const static int max_jet = 1000;
  UInt_t value_jet_n;
  float value_jet_pt[max_jet];
  float value_jet_eta[max_jet];
  float value_jet_phi[max_jet];
  float value_jet_mass[max_jet];
  bool value_jet_looseid[max_jet];
  bool value_jet_tightid[max_jet];
  float value_jet_btag[max_jet];

  // Generator particles
  const static int max_gen = 1000;
  UInt_t value_gen_n;
  float value_gen_pt[max_gen];
  float value_gen_eta[max_gen];
  float value_gen_phi[max_gen];
  float value_gen_mass[max_gen];
  int value_gen_pdgid[max_gen];
  int value_gen_status[max_gen];
};

AOD2NanoAOD::AOD2NanoAOD(const edm::ParameterSet &iConfig)
        : isData(iConfig.getParameter<bool>("isData")) {
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("Events", "Events");

  // Event information
  tree->Branch("run", &value_run);
  tree->Branch("luminosityBlock", &value_lumi_block);
  tree->Branch("event", &value_event);

  // Trigger
  for(size_t i = 0; i < interestingTriggers.size(); i++) {
    tree->Branch(interestingTriggers[i].c_str(), value_trig + i, (interestingTriggers[i] + "/O").c_str());
  }

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
  tree->Branch("Muon_pfRelIso03_all", value_mu_pfreliso03all, "Muon_pfRelIso03_all[nMuon]/F");
  tree->Branch("Muon_pfRelIso04_all", value_mu_pfreliso04all, "Muon_pfRelIso04_all[nMuon]/F");
  tree->Branch("Muon_tightId", value_mu_tightid, "Muon_tightId[nMuon]/O");
  tree->Branch("Muon_softId", value_mu_softid, "Muon_softId[nMuon]/O");
  tree->Branch("Muon_dxy", value_mu_dxy, "Muon_dxy[nMuon]/F");
  tree->Branch("Muon_dxyErr", value_mu_dxyErr, "Muon_dxyErr[nMuon]/F");
  tree->Branch("Muon_dz", value_mu_dz, "Muon_dz[nMuon]/F");
  tree->Branch("Muon_dzErr", value_mu_dzErr, "Muon_dzErr[nMuon]/F");
  tree->Branch("Muon_jetIdx", value_mu_jetidx, "Muon_jetIdx[nMuon]/I");
  tree->Branch("Muon_genPartIdx", value_mu_genpartidx, "Muon_genPartIdx[nMuon]/I");

  // Electrons
  tree->Branch("nElectron", &value_el_n, "nElectron/i");
  tree->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_mass", value_el_mass, "Electron_mass[nElectron]/F");
  tree->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");
  tree->Branch("Electron_pfRelIso03_all", value_el_pfreliso03all, "Electron_pfRelIso03_all[nElectron]/F");
  tree->Branch("Electron_dxy", value_el_dxy, "Electron_dxy[nElectron]/F");
  tree->Branch("Electron_dxyErr", value_el_dxyErr, "Electron_dxyErr[nElectron]/F");
  tree->Branch("Electron_dz", value_el_dz, "Electron_dz[nElectron]/F");
  tree->Branch("Electron_dzErr", value_el_dzErr, "Electron_dzErr[nElectron]/F");
  tree->Branch("Electron_cutBasedId", value_el_cutbasedid, "Electron_cutBasedId[nElectron]/O");
  tree->Branch("Electron_pfId", value_el_pfid, "Electron_pfId[nElectron]/O");
  tree->Branch("Electron_jetIdx", value_el_jetidx, "Electron_jetIdx[nElectron]/I");
  tree->Branch("Electron_genPartIdx", value_el_genpartidx, "Electron_genPartIdx[nElectron]/I");

  // Taus
  tree->Branch("nTau", &value_tau_n, "nTau/i");
  tree->Branch("Tau_pt", value_tau_pt, "Tau_pt[nTau]/F");
  tree->Branch("Tau_eta", value_tau_eta, "Tau_eta[nTau]/F");
  tree->Branch("Tau_phi", value_tau_phi, "Tau_phi[nTau]/F");
  tree->Branch("Tau_mass", value_tau_mass, "Tau_mass[nTau]/F");
  tree->Branch("Tau_charge", value_tau_charge, "Tau_charge[nTau]/I");
  tree->Branch("Tau_decayMode", value_tau_decaymode, "Tau_decayMode[nTau]/I");
  tree->Branch("Tau_chargedIso", value_tau_chargediso, "Tau_chargedIso[nTau]/F");
  tree->Branch("Tau_neutralIso", value_tau_neutraliso, "Tau_neutralIso[nTau]/F");
  tree->Branch("Tau_relIso_all", value_tau_reliso_all, "Tau_relIso_all[nTau]/F");
  tree->Branch("Tau_jetIdx", value_tau_jetidx, "Tau_jetIdx[nTau]/I");
  tree->Branch("Tau_genPartIdx", value_tau_genpartidx, "Tau_genPartIdx[nTau]/I");

  // Photons
  tree->Branch("nPhoton", &value_ph_n, "nPhoton/i");
  tree->Branch("Photon_pt", value_ph_pt, "Photon_pt[nPhoton]/F");
  tree->Branch("Photon_eta", value_ph_eta, "Photon_eta[nPhoton]/F");
  tree->Branch("Photon_phi", value_ph_phi, "Photon_phi[nPhoton]/F");
  tree->Branch("Photon_mass", value_ph_mass, "Photon_mass[nPhoton]/F");
  tree->Branch("Photon_charge", value_ph_charge, "Photon_charge[nPhoton]/I");
  tree->Branch("Photon_pfRelIso03_all", value_ph_pfreliso03all, "Photon_pfRelIso03_all[nPhoton]/F");
  tree->Branch("Photon_jetIdx", value_ph_jetidx, "Photon_jetIdx[nPhoton]/I");
  tree->Branch("Photon_genPartIdx", value_tau_genpartidx, "Photon_genPartIdx[nPhoton]/I");

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
  tree->Branch("Jet_looseId", value_jet_looseid, "Jet_looseId[nJet]/O");
  tree->Branch("Jet_tightId", value_jet_tightid, "Jet_tightId[nJet]/O");
  tree->Branch("Jet_btag", value_jet_btag, "Jet_btag[nJet]/F");

  // Generator particles
  if (!isData) {
    tree->Branch("nGenPart", &value_gen_n, "nGenPart/i");
    tree->Branch("GenPart_pt", value_gen_pt, "GenPart_pt[nGenPart]/F");
    tree->Branch("GenPart_eta", value_gen_eta, "GenPart_eta[nGenPart]/F");
    tree->Branch("GenPart_phi", value_gen_phi, "GenPart_phi[nGenPart]/F");
    tree->Branch("GenPart_mass", value_gen_mass, "GenPart_mass[nGenPart]/F");
    tree->Branch("GenPart_pdgId", value_gen_pdgid, "GenPart_pdgId[nGenPart]/I");
    tree->Branch("GenPart_status", value_gen_status, "GenPart_status[nGenPart]/I");
  }
}

AOD2NanoAOD::~AOD2NanoAOD() {}

void AOD2NanoAOD::analyze(const edm::Event &iEvent,
                          const edm::EventSetup &iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  // Event information
  value_run = iEvent.run();
  value_lumi_block = iEvent.luminosityBlock();
  value_event = iEvent.id().event();

  // Trigger results
  Handle<TriggerResults> trigger;
  iEvent.getByLabel(InputTag("TriggerResults", "", "HLT"), trigger);
  auto psetRegistry = edm::pset::Registry::instance();
  auto triggerParams = psetRegistry->getMapped(trigger->parameterSetID());
  TriggerNames triggerNames(*triggerParams);
  TriggerResultsByName triggerByName(&(*trigger), &triggerNames);
  for (size_t i = 0; i < interestingTriggers.size(); i++) {
    value_trig[i] = false;
  }
  const auto names = triggerByName.triggerNames();
  for (size_t i = 0; i < interestingTriggers.size(); i++) {
    const auto name = interestingTriggers[i];
    auto trig = std::find(names.begin(), names.end(), name);
    if (trig != names.end()) {
      const auto status = triggerByName.state(*trig);
      if (status == 1) {
        value_trig[i] = true;
      }
    }
  }

  // Vertex
  Handle<VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  value_ve_n = vertices->size();
  value_ve_x = vertices->begin()->x();
  value_ve_y = vertices->begin()->y();
  value_ve_z = vertices->begin()->z();
  math::XYZPoint pv(vertices->begin()->position());

  // Muons
  Handle<MuonCollection> muons;
  iEvent.getByLabel(InputTag("muons"), muons);

  value_mu_n = 0;
  const float mu_min_pt = 3;
  std::vector<Muon> selectedMuons;
  for (auto it = muons->begin(); it != muons->end(); it++) {
    if (it->pt() > mu_min_pt) {
      selectedMuons.emplace_back(*it);
      value_mu_pt[value_mu_n] = it->pt();
      value_mu_eta[value_mu_n] = it->eta();
      value_mu_phi[value_mu_n] = it->phi();
      value_mu_charge[value_mu_n] = it->charge();
      value_mu_mass[value_mu_n] = it->mass();
      if (it->isPFMuon() && it->isPFIsolationValid()) {
        auto iso03 = it->pfIsolationR03();
        value_mu_pfreliso03all[value_mu_n] =
            (iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt)/it->pt();
        auto iso04 = it->pfIsolationR04();
        value_mu_pfreliso04all[value_mu_n] =
            (iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/it->pt();
      } else {
        value_mu_pfreliso03all[value_mu_n] = -999;
        value_mu_pfreliso04all[value_mu_n] = -999;
      }
      value_mu_tightid[value_mu_n] = muon::isTightMuon(*it, *vertices->begin());
      value_mu_softid[value_mu_n] = muon::isSoftMuon(*it, *vertices->begin());
      auto trk = it->globalTrack();
      if (trk.isNonnull()) {
        value_mu_dxy[value_mu_n] = trk->dxy(pv);
        value_mu_dz[value_mu_n] = trk->dz(pv);
        value_mu_dxyErr[value_mu_n] = trk->d0Error();
        value_mu_dzErr[value_mu_n] = trk->dzError();
      } else {
        value_mu_dxy[value_mu_n] = -999;
        value_mu_dxyErr[value_mu_n] = -999;
        value_mu_dz[value_mu_n] = -999;
        value_mu_dzErr[value_mu_n] = -999;
      }
      value_mu_genpartidx[value_mu_n] = -1;
      value_mu_jetidx[value_mu_n] = -1;
      value_mu_n++;
    }
  }

  // Electrons
  Handle<GsfElectronCollection> electrons;
  iEvent.getByLabel(InputTag("gsfElectrons"), electrons);

  value_el_n = 0;
  const float el_min_pt = 5;
  std::vector<GsfElectron> selectedElectrons;
  for (auto it = electrons->begin(); it != electrons->end(); it++) {
    if (it->pt() > el_min_pt) {
      selectedElectrons.emplace_back(*it);
      value_el_pt[value_el_n] = it->pt();
      value_el_eta[value_el_n] = it->eta();
      value_el_phi[value_el_n] = it->phi();
      value_el_charge[value_el_n] = it->charge();
      value_el_mass[value_el_n] = it->mass();
      value_el_cutbasedid[value_el_n] = it->passingCutBasedPreselection();
      value_el_pfid[value_el_n] = it->passingPflowPreselection();
      if (it->passingPflowPreselection()) {
        auto iso03 = it->pfIsolationVariables();
        value_el_pfreliso03all[value_el_n] =
            (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/it->pt();
      } else {
        value_el_pfreliso03all[value_el_n] = -999;
      }
      auto trk = it->gsfTrack();
      value_el_dxy[value_el_n] = trk->dxy(pv);
      value_el_dz[value_el_n] = trk->dz(pv);
      value_el_dxyErr[value_el_n] = trk->d0Error();
      value_el_dzErr[value_el_n] = trk->dzError();
      value_el_jetidx[value_el_n] = -1;
      value_el_genpartidx[value_el_n] = -1;
      value_el_n++;
    }
  }

  // Taus
  Handle<PFTauCollection> taus;
  iEvent.getByLabel(InputTag("hpsPFTauProducer"), taus);

  const float tau_min_pt = 15;
  value_tau_n = 0;
  std::vector<PFTau> selectedTaus;
  for (auto it = taus->begin(); it != taus->end(); it++) {
    if (it->pt() > tau_min_pt) {
      selectedTaus.emplace_back(*it);
      value_tau_pt[value_tau_n] = it->pt();
      value_tau_eta[value_tau_n] = it->eta();
      value_tau_phi[value_tau_n] = it->phi();
      value_tau_charge[value_tau_n] = it->charge();
      value_tau_mass[value_tau_n] = it->mass();
      value_tau_decaymode[value_tau_n] = it->decayMode();
      value_tau_chargediso[value_tau_n] = it->isolationPFChargedHadrCandsPtSum();
      value_tau_neutraliso[value_tau_n] = it->isolationPFGammaCandsEtSum();
      value_tau_reliso_all[value_tau_n] = (it->isolationPFChargedHadrCandsPtSum() + it->isolationPFGammaCandsEtSum()) / it->pt();
      value_tau_jetidx[value_tau_n] = -1;
      value_tau_genpartidx[value_tau_n] = -1;
      value_tau_n++;
    }
  }

  // Photons
  Handle<PhotonCollection> photons;
  iEvent.getByLabel(InputTag("photons"), photons);

  value_ph_n = 0;
  const float ph_min_pt = 5;
  std::vector<Photon> selectedPhotons;
  for (auto it = photons->begin(); it != photons->end(); it++) {
    if (it->pt() > ph_min_pt) {
      selectedPhotons.emplace_back(*it);
      value_ph_pt[value_ph_n] = it->pt();
      value_ph_eta[value_ph_n] = it->eta();
      value_ph_phi[value_ph_n] = it->phi();
      value_ph_charge[value_ph_n] = it->charge();
      value_ph_mass[value_ph_n] = it->mass();
      value_ph_pfreliso03all[value_ph_n] = it->ecalRecHitSumEtConeDR03() / it->pt();
      value_ph_jetidx[value_ph_n] = -1;
      value_ph_genpartidx[value_ph_n] = -1;
      value_ph_n++;
    }
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
  // Jet ID recommendations:
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_8_TeV_data_a
  // B-tag recommendations:
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
  Handle<PFJetCollection> jets;
  iEvent.getByLabel(InputTag("ak5PFJets"), jets);
  Handle<JetTagCollection> btags;
  iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);

  const float jet_min_pt = 15;
  value_jet_n = 0;
  std::vector<PFJet> selectedJets;
  for (auto it = jets->begin(); it != jets->end(); it++) {
    if (it->pt() > jet_min_pt) {
      selectedJets.emplace_back(*it);
      value_jet_pt[value_jet_n] = it->pt();
      value_jet_eta[value_jet_n] = it->eta();
      value_jet_phi[value_jet_n] = it->phi();
      value_jet_mass[value_jet_n] = it->mass();
      const auto NHF = it->neutralHadronEnergyFraction();
      const auto NEMF = it->neutralEmEnergyFraction();
      const auto CHF = it->chargedHadronEnergyFraction();
      const auto MUF = it->muonEnergyFraction();
      const auto CEMF = it->chargedEmEnergyFraction();
      const auto NumConst = it->chargedMultiplicity()+it->neutralMultiplicity();
      const auto CHM = it->chargedMultiplicity();
      value_jet_looseid[value_jet_n] = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((std::fabs(it->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || std::fabs(it->eta())>2.4);
      value_jet_tightid[value_jet_n] = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((std::fabs(it->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || std::fabs(it->eta())>2.4);
      value_jet_btag[value_jet_n] = btags->operator[](it - jets->begin()).second;
      value_jet_n++;
    }
  }

  // Generator particles
  if (!isData) {
    Handle<GenParticleCollection> gens;
    iEvent.getByLabel(InputTag("genParticles"), gens);

    value_gen_n = 0;
    std::vector<GenParticle> interestingGenMuons;
    std::vector<GenParticle> interestingGenElectrons;
    std::vector<GenParticle> interestingGenPhotons;
    std::vector<GenParticle> interestingGenTaus;
    for (auto it = gens->begin(); it != gens->end(); it++) {
      const auto status = it->status();
      const auto pdgId = std::abs(it->pdgId());
      if (status == 1 && pdgId == 13) { // muon
        interestingGenMuons.emplace_back(*it);
      }
      if (status == 1 && pdgId == 11) { // electron
        interestingGenElectrons.emplace_back(*it);
      }
      if (status == 1 && pdgId == 22) { // photon
        interestingGenPhotons.emplace_back(*it);
      }
      if (status == 2 && pdgId == 15) { // tau
        interestingGenTaus.emplace_back(*it);
      }
    }

    // Match muons with gen particles and jets
    const auto deltaRMax = 0.3;
    for (auto m = selectedMuons.begin(); m != selectedMuons.end(); m++) {
      // Gen particle matching
      for (auto g = interestingGenMuons.begin(); g != interestingGenMuons.end(); g++) {
        if (deltaR(m->p4(), g->p4()) < deltaRMax) {
          value_gen_pt[value_gen_n] = g->pt();
          value_gen_eta[value_gen_n] = g->eta();
          value_gen_phi[value_gen_n] = g->phi();
          value_gen_mass[value_gen_n] = g->mass();
          value_gen_pdgid[value_gen_n] = g->pdgId();
          value_gen_status[value_gen_n] = g->status();
          value_mu_genpartidx[m - selectedMuons.begin()] = value_gen_n;
          value_gen_n++;
          break;
        }
      }
      // Jet matching
      for(auto j = selectedJets.begin(); j != selectedJets.end(); j++) {
        if (deltaR(m->p4(), j->p4()) < deltaRMax) {
          value_mu_jetidx[m - selectedMuons.begin()] = j - selectedJets.begin();
        }
      }
    }

    // Match electrons with gen particles and jets
    for (auto m = selectedElectrons.begin(); m != selectedElectrons.end(); m++) {
      // Gen particle matching
      for (auto g = interestingGenElectrons.begin(); g != interestingGenElectrons.end(); g++) {
        if (deltaR(m->p4(), g->p4()) < deltaRMax) {
          value_gen_pt[value_gen_n] = g->pt();
          value_gen_eta[value_gen_n] = g->eta();
          value_gen_phi[value_gen_n] = g->phi();
          value_gen_mass[value_gen_n] = g->mass();
          value_gen_pdgid[value_gen_n] = g->pdgId();
          value_gen_status[value_gen_n] = g->status();
          value_el_genpartidx[m - selectedElectrons.begin()] = value_gen_n;
          value_gen_n++;
          break;
        }
      }
      // Jet matching
      for(auto j = selectedJets.begin(); j != selectedJets.end(); j++) {
        if (deltaR(m->p4(), j->p4()) < deltaRMax) {
          value_el_jetidx[m - selectedElectrons.begin()] = j - selectedJets.begin();
        }
      }
    }

   // Match photons with gen particles and jets
    for (auto m = selectedPhotons.begin(); m != selectedPhotons.end(); m++) {
      // Gen particle matching
      for (auto g = interestingGenPhotons.begin(); g != interestingGenPhotons.end(); g++) {
        if (deltaR(m->p4(), g->p4()) < deltaRMax) {
          value_gen_pt[value_gen_n] = g->pt();
          value_gen_eta[value_gen_n] = g->eta();
          value_gen_phi[value_gen_n] = g->phi();
          value_gen_mass[value_gen_n] = g->mass();
          value_gen_pdgid[value_gen_n] = g->pdgId();
          value_gen_status[value_gen_n] = g->status();
          value_ph_genpartidx[m - selectedPhotons.begin()] = value_gen_n;
          value_gen_n++;
          break;
        }
      }
      // Jet matching
      for(auto j = selectedJets.begin(); j != selectedJets.end(); j++) {
        if (deltaR(m->p4(), j->p4()) < deltaRMax) {
          value_ph_jetidx[m - selectedPhotons.begin()] = j - selectedJets.begin();
        }
      }
    }

    // Match taus with gen particles and jets
    for (auto m = selectedTaus.begin(); m != selectedTaus.end(); m++) {
      // Gen particle matching
      for (auto g = interestingGenTaus.begin(); g != interestingGenTaus.end(); g++) {
        // Subtract neutrinos from generator particle
        auto p4 = g->p4();
        for (auto d = g->begin(); d != g->end(); d++) {
          const auto pdgId = d->pdgId();
          if (pdgId == 12 || pdgId == 14 || pdgId == 16 || pdgId == 18) {
            p4 -= d->p4();
          }
        }
        if (deltaR(m->p4(), p4) < deltaRMax) {
          value_gen_pt[value_gen_n] = g->pt();
          value_gen_eta[value_gen_n] = g->eta();
          value_gen_phi[value_gen_n] = g->phi();
          value_gen_mass[value_gen_n] = g->mass();
          value_gen_pdgid[value_gen_n] = g->pdgId();
          value_gen_status[value_gen_n] = g->status();
          value_tau_genpartidx[m - selectedTaus.begin()] = value_gen_n;
          value_gen_n++;
          break;
        }
      }
      // Jet matching
      for(auto j = selectedJets.begin(); j != selectedJets.end(); j++) {
        if (deltaR(m->p4(), j->p4()) < deltaRMax) {
          value_tau_jetidx[m - selectedTaus.begin()] = j - selectedJets.begin();
        }
      }
    }

  } // !isData

  // Fill event
  tree->Fill();
}

void AOD2NanoAOD::beginJob() {}

void AOD2NanoAOD::endJob() {}

DEFINE_FWK_MODULE(AOD2NanoAOD);

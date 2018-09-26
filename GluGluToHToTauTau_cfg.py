import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("AOD2NanoAOD")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.categories.append("AOD2NanoAOD")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

# Set the maximum number of events to be processed (-1 processes all events)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

# Define files of dataset
files = FileUtils.loadListFromFile(
    "data/CMS_MonteCarlo2012_Summer12_DR53X_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6-tauPolarOff_AODSIM_PU_S10_START53_V19-v1_00000_file_index.txt")
process.source = cms.Source(
    "PoolSource", fileNames=cms.untracked.vstring(*files))

# Number of events to be skipped (0 by default)
process.source.skipEvents = cms.untracked.uint32(0)

# Register fileservice for output file
process.aod2nanoaod = cms.EDAnalyzer("AOD2NanoAOD")
process.TFileService = cms.Service(
    "TFileService", fileName=cms.string("GluGluToHToTauTau.root"))

process.p = cms.Path(process.aod2nanoaod)

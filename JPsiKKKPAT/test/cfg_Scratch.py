import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False)
)
runOnMC=False
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200))

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
'/store/data/Run2017B/Charmonium/AOD/PromptReco-v1/000/297/050/00000/AC4BE68E-3F56-E711-9892-02163E01A1B7.root'
))

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4', '')

#select generaltracks
process.oniaSelectedTracks=cms.EDFilter("TrackSelector", 
      src = cms.InputTag("generalTracks"), 
      cut = cms.string('pt > 0.7 && abs(eta) <= 3.0' 
                                '&& charge !=0' 
                                '&& quality(\"highPurity\")'
                                ) 
)

#convert to PAT
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.CandidateSelectedTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
            src=cms.InputTag("oniaSelectedTracks"),
            particleType=cms.string('K+')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

#make muon with trigger matching embedded
process.load('HeavyFlavorAnalysis.Onia2MuMu.oniaPATMuonsWithTrigger_cff')
#select them, adjust as needed, currently soft muon selection with pT>3.5 and |eta|<1.9
process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('oniaPATMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
   ),
   filter = cms.bool(True)
)

process.mkcands = cms.EDAnalyzer("JPsiKKKPATRunII",
     HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
     beamSpotTag = cms.InputTag("offlineBeamSpot"),
     VtxSample   = cms.InputTag('offlinePrimaryVertices'),
     Trak        = cms.InputTag('patSelectedTracks'), #selectedPatTracks
     muons       = cms.InputTag('oniaSelectedMuons'), #oniaPATMuonsWithoutTrigger-  #selectedPatMuons
     revtxtrks   = cms.InputTag('generalTracks'),
     trackQualities = cms.untracked.vstring('loose','tight','highPurity'),
     TriggersForMatching = cms.untracked.vstring("HLT_Dimuon16_Jpsi_v3", "HLT_DoubleMu4_JpsiTrk_Displaced_v3", "HLT_Dimuon6_Jpsi_NoVertexing_v3", "HLT_DoubleMu4_3_Jpsi_Displaced_v3"), #HLT_DoubleMu4_3_Jpsi_Displaced_v2
     FiltersForMatching = cms.untracked.vstring("hltDisplacedmumuFilterDimuon16Jpsi", "hltJpsiTkVertexFilter", "hltDimuon6JpsiL3Filtered", "hltDisplacedmumuFilterDoubleMu43Jpsi" ) #hltDisplacedmumuFilterDoubleMu43Jpsi
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuOniaRun2012_PhaseSpMC2016_ntpl.root')
)
process.mySequence = cms.Sequence(
                   process.oniaPATMuonsWithTriggerSequence *
                   process.oniaSelectedMuons *
                   process.oniaSelectedTracks *
                   process.CandidateSelectedTracks *
                   process.patSelectedTracks *
                   process.mkcands
)
process.p = cms.Path(process.mySequence)
process.schedule = cms.Schedule(process.p)

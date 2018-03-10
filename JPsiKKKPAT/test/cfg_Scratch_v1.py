import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False)
)
MC=False


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


process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
                            '/store/data/Run2017B/Charmonium/AOD/PromptReco-v1/000/297/050/00000/AC4BE68E-3F56-E711-9892-02163E01A1B7.root'
                            ))


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4', '')

process.load("PhysicsTools.PatAlgos.patSequences_cff")


############################
### Add track candidates ###
############################
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
# Require NOT to check overlap with muons and electrons
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)


#########################
### Set up PAT tracks ###
#########################
from PhysicsTools.PatAlgos.tools.trackTools import *
from PhysicsTools.PatAlgos.tools.coreTools  import *

makeTrackCandidates(process,                  # patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
        particleType='K+',                   # particle type (for assigning a mass), changed by yik to K+ from pi+
        preselection='pt > 0.5',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available, changed to 0.5 from 0.1 by yik
        selection='pt > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands'), changed to 0.5 from 0.1 by yik
        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
        isoDeposits=[],
        mcAs=None                           # Replicate MC match as the one used for Muons
        );                                    # you can specify more than one collection for this
l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

removeMCMatching(process, ['All'], outputModules = [])


### Trigger matching ###

process.cleanMuonHLTriggerMatch = cms.EDProducer(
		'PATTriggerMatcherDRLessByR',
		src                   = cms.InputTag('cleanPatMuons'),
		matched               = cms.InputTag('patTrigger'),
		matchedCuts           = cms.string('path("HLT_DoubleMu4_3_Jpsi_Displaced_v*") || path("HLT_DoubleMu4_JpsiTrk_Displaced_v*")'),
		maxDeltaR             = cms.double(0.1),
		resolveAmbiguities    = cms.bool(True),
		resolveByMatchQuality = cms.bool(True))


### Switch on PAT trigger ###
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonHLTriggerMatch'], hltProcess = triggerProcessName, outputModule = '')


# turn off MC matching for the process
process.patDefaultSequence.remove(process.muonMatch)
process.patDefaultSequence.remove(process.electronMatch)
process.patDefaultSequence.remove(process.photonMatch)
process.patDefaultSequence.remove(process.tauMatch)
process.patDefaultSequence.remove(process.tauGenJetMatch)

process.patDefaultSequence.remove(process.tauGenJets)
process.patDefaultSequence.remove(process.tauGenJetsSelectorAllHadrons)

process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonsLegacy)            
process.patDefaultSequence.remove(process.patJetPartonAssociationLegacy)  
process.patDefaultSequence.remove(process.patJetFlavourAssociationLegacy) 
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)

process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.patDefaultSequence.remove(process.selectedPatTaus)
process.patDefaultSequence.remove(process.cleanPatTaus)
process.patDefaultSequence.remove(process.countPatTaus)

process.patDefaultSequence.remove(process.selectedPatElectrons)
process.patDefaultSequence.remove(process.cleanPatElectrons)
process.patDefaultSequence.remove(process.countPatElectrons)

process.patDefaultSequence.remove(process.selectedPatPhotons)
process.patDefaultSequence.remove(process.cleanPatPhotons)
process.patDefaultSequence.remove(process.countPatPhotons)

process.patDefaultSequence.remove(process.patMETs)

##################################good collisions############################################
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                                      minimumNDOF = cms.uint32(4) ,
                                                      maxAbsZ = cms.double(24),
                                                      maxd0 = cms.double(2)
)

process.noscraping = cms.EDFilter("FilterOutScraping",
applyfilter = cms.untracked.bool(True),
debugOn = cms.untracked.bool(False),
numtrack = cms.untracked.uint32(10),
thresh = cms.untracked.double(0.25)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Prune generated particles to muons and their parents
process.genMuons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *  ",                     # this is the default
            "++keep abs(pdgId) = 13",        # keep muons and their parents
            "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
      )
 )
 
process.mkcands = cms.EDAnalyzer('JPsiKKKPATRunII',
     HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
     beamSpotTag = cms.InputTag("offlineBeamSpot"),
     primaryVertexTag   = cms.InputTag('offlinePrimaryVertices'),
     DoMonteCarloTree = cms.untracked.bool(MC), # change the MC bool tag in the beginning correspondingly!
     muons       = cms.InputTag('cleanPatMuonsTriggerMatch'), #selectedPatMuons
     Trak        = cms.InputTag('cleanPatTrackCands'), 
     revtxtrks   = cms.InputTag('generalTracks'),
     inputGEN = cms.InputTag("genParticles"),
     TriggersForMatching = cms.untracked.vstring( "HLT_DoubleMu4_3_Jpsi_Displaced_v2", "HLT_DoubleMu4_3_Jpsi_Displaced_v3", "HLT_DoubleMu4_3_Jpsi_Displaced_v4", "HLT_DoubleMu4_3_Jpsi_Displaced_v5", "HLT_DoubleMu4_3_Jpsi_Displaced_v7","HLT_DoubleMu4_JpsiTrk_Displaced_v2",  "HLT_DoubleMu4_JpsiTrk_Displaced_v3", "HLT_DoubleMu4_JpsiTrk_Displaced_v4", "HLT_DoubleMu4_JpsiTrk_Displaced_v5", "HLT_DoubleMu4_JpsiTrk_Displaced_v7"),
     FiltersForMatching = cms.untracked.vstring("hltDisplacedmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi" ), #hltDisplacedmumuFilterDoubleMu43Jpsi
#     resolvePileUpAmbiguity = cms.untracked.bool(True),
#     addBlessPrimaryVertex = cms.untracked.bool(True),
#     DoJPsiMassConstraint = cms.untracked.bool(False),
#     UseBpKDR = cms.untracked.bool(False),
#     UseJpsiKDR = cms.untracked.bool(True)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuOnia_Run2_JPsiKKKPAT_ntpl.root')
)                                   
                                
### INITIALIZE ###
process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

if( not MC ):
  process.mySequence = cms.Sequence(
                   process.filter *
                   process.patDefaultSequence *
                   process.mkcands
)
else:
        process.mySequence = cms.Sequence(
                   process.mkcands
)

process.p = cms.Path(process.mySequence)
process.schedule = cms.Schedule(process.p)

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4') #2022 C D E F 
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v10') #2022 G
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_HLT_v4')
#process.GlobalTag = GlobalTag(process.GlobalTag, '120X_mcRun3_2021_realistic_v5')
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_postEE_v1')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

from mcList import inputList

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      inputList[1500:2000]
#'file:/asanchez/data/store/data/Run2016G/Charmonium/MINIAOD/23Sep2016-v1/A4B4AC67-B996-E611-9ECD-008CFAFBE8CE.root',
#MiniAOD
#'/store/data/Run2018C/Charmonium/MINIAOD/PromptReco-v2/000/319/756/00000/EEF6CEC1-698B-E811-8081-02163E00AF5F.root',
#'/store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/355/872/00000/18146588-e660-4b01-80b6-7611087b65a7.root',
#'/store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/356/076/00000/286c5c07-05d2-4cec-acc0-05514991700a.root',
#'/store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/357/482/00000/df4fd089-0144-4bdd-9613-7b6cc1009841.root',
#'/store/data/Run2022F/ParkingDoubleMuonLowMass1/MINIAOD/PromptReco-v1/000/360/390/00000/0bc48cdc-66b8-4316-84ae-f2934a2304fb.root',
#'/store/data/Run2022B/DoubleMuonLowMass/MINIAOD/PromptReco-v1/000/355/558/00000/895a4444-ac75-4ab6-a682-b2580508dc5d.root',       
#'/store/mc/Run3Summer21MiniAOD/BcToJpsiPi_TuneCP5_14TeV_pythia8_evtgen/MINIAODSIM/Pilot_120X_mcRun3_2021_realistic_v5_ext1-v2/70000/3b915027-ecc8-40ab-b589-dca511dc7e8b.root'
 )
)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_DoubleMu4_JpsiTrk_Bc_v*',),
								#	'HLT_DoubleMu4_MuMuTrk_Displaced_v*',                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.load("myAnalyzers.JPsiKsPAT.slimmedMuonsTriggerMatcher_cfi")  

process.load("myAnalyzers.JPsiKsPAT.JPsiPiRootupler_cfi")
#process.rootuple.dimuons = cms.InputTag('slimmedMuonsWithTrigger') 

process.TFileService = cms.Service("TFileService",

       fileName = cms.string('Rootuple_Bc-MiniAOD-MC-triggerOn-4.root'),
  
)

process.mySequence = cms.Sequence(
                                   process.triggerSelection *
#     				   process.slimmedMuonsWithTriggerSequence *
                                   process.rootuple
				   )

process.p = cms.Path(process.mySequence)

#process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
#process.p = cms.Path(process.rootuple)



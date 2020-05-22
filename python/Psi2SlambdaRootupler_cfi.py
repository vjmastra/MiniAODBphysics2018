import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('Psi2SLambda',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secundaryVerticesPtr = cms.InputTag("slimmedLambdaVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerInput = cms.InputTag("slimmedPatTrigger"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(True),
                          OnlyGen = cms.bool(False),
                          )

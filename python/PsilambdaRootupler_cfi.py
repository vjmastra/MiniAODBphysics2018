import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiLambda',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secundaryVerticesPtr = cms.InputTag("slimmedLambdaVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#                          unpackFilterLabels = cms.bool(True),
                          objects = cms.InputTag("unpackedPatTrigger"),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(True),
                          OnlyGen = cms.bool(False),
                          )

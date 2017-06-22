import FWCore.ParameterSet.Config as cms

flashggTorchDumperBarrel = cms.EDAnalyzer('FlashggPhoIdDumperBarrel',
                                          diphotonsInput = cms.InputTag('flashggPreselectedDiPhotons'),
                                          # diphotonsInput = cms.InputTag('flashggUpdatedIdMVADiPhotons'),
                                          
                                          output = cms.untracked.string("rechits-barrel.t7"),
                                          windowHalfWidth = cms.untracked.uint32(17),
                                          windowHalfHeight = cms.untracked.uint32(17),
                                          writeSparse = cms.untracked.bool(True),

                                          photonIdInputVarsInputTag = cms.InputTag("flashggUpdatedIdMVADiPhotons"),
                                          writePhotonIdInputVars = cms.untracked.bool(True),

                                          normalizeRecHitsToMax = cms.untracked.bool(False),

                                          vertexCandidateMapTag = cms.InputTag("flashggVertexMapNonUnique"),
                                          )

flashggTorchDumperEndcap = cms.EDAnalyzer('FlashggPhoIdDumperEndcap',
                                          diphotonsInput = cms.InputTag('flashggPreselectedDiPhotons'),
                                          # diphotonsInput = cms.InputTag('flashggUpdatedIdMVADiPhotons'),
                                          
                                          output = cms.untracked.string("rechits-endcap.t7"),
                                          windowHalfWidth = cms.untracked.uint32(17),
                                          windowHalfHeight = cms.untracked.uint32(17),
                                          writeSparse = cms.untracked.bool(True),

                                          photonIdInputVarsInputTag = cms.InputTag("flashggUpdatedIdMVADiPhotons"),
                                          writePhotonIdInputVars = cms.untracked.bool(True),

                                          normalizeRecHitsToMax = cms.untracked.bool(False),

                                          vertexCandidateMapTag = cms.InputTag("flashggVertexMapNonUnique"),
                                          )

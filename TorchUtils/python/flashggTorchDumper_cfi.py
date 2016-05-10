import FWCore.ParameterSet.Config as cms

flashggTorchDumper = cms.EDAnalyzer('FlashggTorchDumper',
                                    diphotonsInput = cms.InputTag('flashggPreselectedDiPhotons'),
                                    # diphotonsInput = cms.InputTag('flashggUpdatedIdMVADiPhotons'),

                                    barrel = cms.untracked.PSet(
                                        output = cms.untracked.string("rechits-barrel.t7"),
                                        windowHalfWidth = cms.untracked.uint32(3),
                                        windowHalfHeight = cms.untracked.uint32(11),
                                        writeSparse = cms.untracked.bool(True),
                                        ),

                                    endcap = cms.untracked.PSet(
                                        output = cms.untracked.string("rechits-endcap.t7"),
                                        windowHalfWidth = cms.untracked.uint32(3),
                                        windowHalfHeight = cms.untracked.uint32(11),
                                        writeSparse = cms.untracked.bool(True),
                                        ),
                                   )

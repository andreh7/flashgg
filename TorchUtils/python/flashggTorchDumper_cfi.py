import FWCore.ParameterSet.Config as cms

flashggTorchDumper = cms.EDAnalyzer('FlashggTorchDumper',
                                    diphotonsInput = cms.InputTag('flashggPreselectedDiPhotons'),
                                    # diphotonsInput = cms.InputTag('flashggUpdatedIdMVADiPhotons'),
                                    barrelOutput = cms.untracked.string("rechits-barrel.t7"),
                                    endcapOutput = cms.untracked.string("rechits-endcap.t7"),
                                  )

import FWCore.ParameterSet.Config as cms

flashggTorchDumper = cms.EDAnalyzer('FlashggTorchDumper',
                                    diphotonsInput = cms.InputTag('flashggPreselectedDiPhotons'),
                                  )

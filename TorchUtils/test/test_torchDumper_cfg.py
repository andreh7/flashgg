#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import os

# SYSTEMATICS SECTION
dropVBFInNonGold = True

process = cms.Process("FLASHggTorchDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
else:
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)
if dropVBFInNonGold:
    process.flashggVBFTag.SetArbitraryNonGoldMC = True
    process.flashggVBFTag.DropNonGoldData = True
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

systlabels = [""]
phosystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels = []

# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()
print "customize.processId:",customize.processId
# load appropriate scale and smearing bins here
# systematics customization scripts will take care of adjusting flashggDiPhotonSystematics
#process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")

# Or use the official tool instead
useEGMTools(process)

customizeSystematicsForBackground(process)


# cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,jetsystlabels,jetSystematicsInputTags)

###### Dumper section

from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
        #                             "file:myMicroAODOutputFile.root"
        # "/store/user/holzner/flashgg/RunIIFall15DR76-1_3_0-25ns_ah/1_3_0/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/RunIIFall15DR76-1_3_0-25ns_ah-1_3_0-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160502_144226/0000/myMicroAODOutputFile_10.root",

        # GJet
        '/store/user/holzner/flashgg/RunIIFall15DR76-1_3_0-25ns_ah/1_3_0/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIIFall15DR76-1_3_0-25ns_ah-1_3_0-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160502_160042/0000/myMicroAODOutputFile_100.root',

        # "root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/VBFHToGG_M-120_13TeV_powheg_pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160210_045711/0000/myMicroAODOutputFile_1.root"
))

#----------
# rechits dumper
#----------
if True:
    process.load("flashgg.TorchUtils.flashggTorchDumper_cfi") 

    from flashgg.MetaData.JobConfig import customize
    customize.parse()

    if customize.options.has_key('jobId'):

        assert hasattr(customize, 'processId')
        assert customize.processId != ""
        import re

        for subdet in [ process.flashggTorchDumperBarrel,
                        process.flashggTorchDumperEndcap ]:
            subdet.output = cms.untracked.string(str(customize.processId) + "_" + re.sub("\.t7$","", subdet.output.value()) + "_%d.t7" % customize.options.jobId)
    else:
        raise Exception("jobId not found")


#----------

# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v1",
                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter

# Split WH and ZH
process.genFilter = cms.Sequence()
if (customize.processId.count("wh") or customize.processId.count("zh")) and not customize.processId.count("wzh"):
    process.load("flashgg/Systematics/VHFilter_cfi")
    process.genFilter += process.VHFilter
    process.VHFilter.chooseW = bool(customize.processId.count("wh"))
    process.VHFilter.chooseZ = bool(customize.processId.count("zh"))

# Split out prompt-fake or fake-fake
process.finalFilter = cms.Sequence()
if (customize.processId.count("qcd") or customize.processId.count("gjet")) and customize.processId.count("fake"):
    process.load("flashgg/Systematics/PromptFakeFilter_cfi")
    process.finalFilter += process.PromptFakeFilter
    if (customize.processId.count("promptfake")):
        process.PromptFakeFilter.doPromptFake = cms.bool(True)
        process.PromptFakeFilter.doFakeFake =cms.bool(False)
    elif (customize.processId.count("fakefake")):
        process.PromptFakeFilter.doPromptFake =cms.bool(False)
        process.PromptFakeFilter.doFakeFake =cms.bool(True)
    else:
        raise Exception,"Mis-configuration of python for prompt-fake filter"

process.p = cms.Path(process.dataRequirements*
                     process.genFilter*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggMuonSystematics*process.flashggElectronSystematics*
                     (process.flashggUnpackedJets*process.jetSystematicsSequence)*
                     (process.flashggTagSequence*process.systematicsTagSequences)*
                     process.flashggSystTagMerger*
                     process.finalFilter *
                     process.flashggTorchDumperBarrel * 
                     process.flashggTorchDumperEndcap
                     )

print "--- Dumping modules that take diphotons as input: ---"
mns = process.p.moduleNames()
for mn in mns:
    module = getattr(process,mn)
    if hasattr(module,"src") and type(module.src) == type(cms.InputTag("")) and module.src.value().count("DiPhoton"):
        print str(module),module.src
    elif hasattr(module,"DiPhotonTag"):
        print str(module),module.DiPhotonTag
print
printSystematicInfo(process)


##############
## Dump EDM ##
##############

# process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('CustomizeWillChangeThisAnyway.root'),
#                                outputCommands = cms.untracked.vstring('keep *') # dump everything! small tests only!
#                                )
# process.e = cms.EndPath(process.out)


# set default options if needed
customize.setDefault("maxEvents",300)
customize.setDefault("targetLumi",2.61e+3)
# call the customization
customize(process)

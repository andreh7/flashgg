#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import os

execfile(os.path.expandvars("$CMSSW_BASE/src/flashgg/Systematics/test/workspaceStd.py"))

#----------------------------------------------------------------------
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
            subdet.output = cms.untracked.string(str(customize.processId) + "_" + re.sub("\.t7$","", subdet.output.value()) + "_%d.npz" % customize.options.jobId)
    else:
        raise Exception("jobId not found")

#----------
# add another instance of the diphoton update producer 
#----------
# take the standard one first
process.load("flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi")

### process.flashggPreselectedDiPhotons2 = cms.EDProducer("FlashggDiPhotonWithUpdatedPhoIdMVAProducer",
###                                                src                      = cms.InputTag("flashggPreselectedDiPhotons"),
###                                                rhoFixedGridCollection   = cms.InputTag('fixedGridRhoAll'),
###                                                #photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_r9shift_barrel.xml"),
###                                                #photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_barrel.xml"),
###                                                photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_r9s4EtaWshift_barrel.xml"),
###                                                photonIdMVAweightfile_EE = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_76X_25ns_endcap.xml"),
###                                                # comment this out to disable all corrections in this module
###                                                # correctionFile           = cms.FileInPath("flashgg/MicroAOD/data/transformation_76X_v2.root"),
###                                                Debug                    = cms.bool(False),
### 
###                                                # do not correct a second time
###                                                applyCorrections         = cms.bool(False),
###                                               )

process.flashggPreselectedDiPhotons2 = process.flashggUpdatedIdMVADiPhotons.clone()


# do not correct a second time, disable corrections
process.flashggPreselectedDiPhotons2.applyCorrections = cms.bool(False)
# process.flashggPreselectedDiPhotons2.correctionFile = cms.FileInPath("")
# process.flashggPreselectedDiPhotons2.correctionFile = cms.FileInPath("/dev/null")
del process.flashggPreselectedDiPhotons2.correctionFile

# use the first one as source. Note that this module has a similar 
# name but the type is different
process.flashggPreselectedDiPhotons2.src = cms.InputTag("flashggPreselectedDiPhotons")

#----------
# add our own dumpers after finalFilter
#----------
pos = process.p.index(process.finalFilter)
process.p.insert(pos + 1, process.flashggPreselectedDiPhotons2)
process.p.insert(pos + 2, process.flashggTorchDumperBarrel)
process.p.insert(pos + 3, process.flashggTorchDumperEndcap)

# remove the actual tagsDumper
process.p.remove(process.tagsDumper)

# take the collection produced based on the preselected photons (without applying corrections again)
# instead of the collection produced just with the corrections (but without preselection)
for mod in (process.flashggTorchDumperBarrel, process.flashggTorchDumperEndcap):

    mod.diphotonsInput            = cms.InputTag('flashggPreselectedDiPhotons2')
    mod.photonIdInputVarsInputTag = cms.InputTag('flashggPreselectedDiPhotons2')

    # mod.diphotonsInput = cms.InputTag('flashggUpdatedIdMVADiPhotons')

    # mod.writePhotonIdInputVars = cms.untracked.bool(False)


#----------
# from Inna
#----------
process.flashggPreselectedDiPhotons.cut=cms.string(
        "mass > 95"
        " && leadingPhoton.pt > 18 && subLeadingPhoton.pt > 18"
        " && abs(leadingPhoton.superCluster.eta)<2.5 && abs(subLeadingPhoton.superCluster.eta)<2.5 "
        " && ( abs(leadingPhoton.superCluster.eta)<1.4442 || abs(leadingPhoton.superCluster.eta)>1.566)"
        " && ( abs(subLeadingPhoton.superCluster.eta)<1.4442 || abs(subLeadingPhoton.superCluster.eta)>1.566)"
        " && (leadingPhoton.pt > 14 && leadingPhoton.hadTowOverEm < 0.15 && (leadingPhoton.full5x5_r9>0.8 || leadingPhoton.chargedHadronIso<20 || leadingPhoton.chargedHadronIso<(0.3*leadingPhoton.pt)))"
        " && (subLeadingPhoton.pt > 14 && subLeadingPhoton.hadTowOverEm < 0.15 && (subLeadingPhoton.full5x5_r9>0.8 || subLeadingPhoton.chargedHadronIso<20 || subLeadingPhoton.chargedHadronIso<(0.3*subLeadingPhoton.pt)))"
        )
#----------

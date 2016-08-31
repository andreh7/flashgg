import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

import os

execfile(os.path.expandvars("${CMSSW_BASE}/src/flashgg/MicroAOD/test/microAODstd.py"))

process.out.outputCommands.append('drop patPackedCandidates_packedPFCandidates_*_*')

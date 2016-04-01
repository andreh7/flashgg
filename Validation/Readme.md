MICROAOD Production
================
```
cd MetaData/work
```
Prepare a json file like the following (myjson.json):
```
{
    "data" : [],
    "sig"  : [],
    "bkg"  : ["/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM",
              "/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM"]
}
```
Run the following command (-h or --help option available):
```
prepareCrabJobs.py -p ../../MicroAOD/test/microAODstd.py -s gjet.json --mkPilot -C A_NAME_FOR_THE_CAMPAIGN -O T2_US_UCSD -o /store/user/sani/flashgg -L 10 --lumiMask https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON.txt -U 2
```
Before running the jobs take a look at the crabConfig_TEMPLATE.py file to see if the options are ok with you.
```
voms-proxy-init -voms cms
cd A_NAME_FOR_THE_CAMPAIGN
echo crabConfig_*.py | xargs -n 1 crab sub
```
Once you have the microAOD you have to produce the catalogue corresponding to your production.
This is MANDATORY to run on the MICROAOD with flashgg !
To make the catalogue please refer to the instructions available here:
https://github.com/cms-analysis/flashgg/tree/master/MetaData (README file)

Processing MICROAOD
===================
- opttree production

You have to edit the json file to specify which samples you are running on, something like:
```
{
    "processes" : {
        "signal" : [
        "/GluGluHToGG_M-125_13TeV_powheg_pythia8",
        "/GluGluHToGG_M-120_13TeV_powheg_pythia8"
        ],
        "data" : [ 
        "/DoubleEG/sethzenz-RunIISpring15-Prompt-1_1_0-25ns-1_1_0-v0-Run2015D-PromptReco-v4-7b05e4f729d8444b2665263390b46268/USER"
        ]
     },
     "cmdLine" : "campaign=RunIIFall15DR76-1_3_0-25ns_ext1 targetLumi=2.5e+3 puTarget=1.34e+05,6.34e+05,8.42e+05,1.23e+06,2.01e+06,4.24e+06,1.26e+07,4.88e+07,1.56e+08,3.07e+08,4.17e+08,4.48e+08,4.04e+08,3.05e+08,1.89e+08,9.64e+07,4.19e+07,1.71e+07,7.85e+06,4.2e+06,2.18e+06,9.43e+05,3.22e+05,8.9e+04,2.16e+04,5.43e+03,1.6e+03,551,206,80.1,31.2,11.9,4.38,1.54,0.518,0.165,0.0501,0.0144,0.00394,0.00102,0.000251,5.87e-05,1.3e-05,2.74e-06,5.47e-07,1.04e-07,1.86e-08,3.18e-09,5.16e-10,9.35e-11"
}
```
  - IMPORTANT parameters:
    - campaing, it refers to a particular version of the MICROAOD, you need to have the corresponding catalogue in MetaData/data/campaing_name/dataset.json
    - puTarget, only when running on MC. The official array will be shared somehow among the group (to get your own please look below)
    - Beware that the dataset names have to match the ones in MetaData/data/cross_sections.json and in the catalogue.

Then you can run the following command:
```
fggRunJobs.py --load your_own_json.json -H -D -P -n how_many_jobs_you_want -d local_output_dir -x cmsRun opttreeDumper.py maxEvents=-1 -q batch_queue --no-use-tarball
```
If everything goes ok you will happily see your trees...

Getting puTarget
================
```
pileupCalc.py lumi -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260426_13TeV_PromptReco_Collisions15_25ns_JSON.txt --inputLumiJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt  --calcMode=“true" --minBiasXsec 72000 --maxPileupBin 52 --numPileupBins 52 myHist.root
```
(with updated parameters). Then you can dump the bin content of myHist.root (no particular normalization is needed).

ZeeValidation
=============
The first step is like the previous one, just replace the json with the proper one (e.g. zee_single.json) and the configuration to be run to zeeValidationDumper.py.
Once you have the corresponding data and MC trees you can produce plots following the instructions below (no CMSSW is needed here).
```
git clone https://github.com/matteosan1/plotZeeValidation
cd plotZeeValidation
root -l 
.L plotter.C++
plotter(data_file_name, mc_file_name)
python combinePlots.py output_file_of_previous_step
```
This should create a .root file with all the plots.
The plots (and the event categorization) are defined in the following files:
- plotvariables.dat (histogram definition)
- categories.dat (category definition)
- selection.dat (event selection definition as a ROOT TCut string)

The final python script combines data and MC histos and put them together with the appropriate style.

Deriving Shower Shape Corrections
=================================
In order to derive shower shape transformation you can use makeTransformation.py script in https://github.com/matteosan1/plotZeeValidation.
Once the Zee validation ntuples have been produced you can run the script three times with different parameters to: 
1) produce the .root file with the necessary plots (python makeTransformation.py -p mc_ntuple_name.root target_ntuple_name.root)
2) derive the transformations (python makeTransformation.py -c)
3) test the transformations on the original ntuples (python makeTransformation.py -t mc_ntuple_name.root target_ntuple_name.root)

IDMVA with systematic plot
==========================
The necessary histograms used as inputs are produced in the ZeeValidation step. More generally you need to add there all the 
plots for nominal, up, down variations which you want to plot.
At this point you need to modify the david_template.py file in order to able to read the name of these plots.
Then:
```
python morphbands,py david_template.py yourROOTFile.root
```
The output is another .root file with plots nomeplot_bottom nomeplot_top of the corrected histograms.
Finally with the idmva_syst.py macro the systematic band is added to the distribution.

In order to apply the additional linearly increasing systematic you can use the linearCorrection.py script (available in the same directory). Beware you need to run it BEFORE the morphbands step described above.
To run it just choose the parameter of the linear correction (xmin, xmax, ymin, ymax) and update the input filename.

TnP Measurement within flashgg
==============================
BEWARE: NEED TO UPDATE THE GEN CONTENT IN FLASHGG TO AVOID USING STATUS 23 ELECTRONS !!!

EGM TnP tool is fully integrated in flashgg. For more details about it please look at https://twiki.cern.ch/twiki/bin/view/CMSPublic/ElectronTagAndProbe
The following steps are going to describe how to produce the necessary trees and how to fit the corresponding invariant mass distributions.

- Preselection
From flashgg/Validation/test:
```
fggRunJobs.py --load tnpMC.json -n 50 -d ~/work/tnpmc -x cmsRun makeTreePhotonsMC.py maxEvents=-1 --no-use-tarball -q 1nd
```
(before running please check that in ../python/treeMakerOptionsPhotons_cfi.py the id_cut parameter of process.goodPhotonTagsIDMVA (0.6) and process.goodPhotonProbesIDMVA (-0.9) are correct) Then run a similar command for data.
Once the jobs are completed:
```
hadd  TnPTree_mc_newID.root ~/work/tnpmc/*.root 
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/getTemplatesFromMC.py --input TnPTree_mc_newID.root -o templatePreselEB.root -d PhotonToRECO --idprobe passingPresel --var1Name probe_Pho_r9 --var1Bins 0,0.85,1 --var2Name probe_sc_abseta --var2Bins 0,1.479
(write down MC efficiencies printed out by the script)
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/getTemplatesFromMC.py --input TnPTree_mc_newID.root -o templatePreselEE.root -d PhotonToRECO --idprobe passingPresel --var1Name probe_Pho_r9 --var1Bins 0,0.90,1 --var2Name probe_sc_abseta --var2Bins 1.479,2.5
(write down MC efficiencies printed out by the script)
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitPreselEB.py --idLabel passingPresel -t templatePreselEB.root --var1Bins 0,0.85,1 --var2Bins 0,1.479
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitPreselEE.py --idLabel passingPresel -t templatePreselEE.root --var1Bins 0,0.90,1 --var2Bins 1.479,2.5
```
Before running the following commands check the binning defined in the fitterWithManyTemplatesPresel.py configuration (and include the appropriate include for EB or EE).
```
cmsRun fitterWithManyTemplatesPresel.py isMC=False inputFileName=TnPTree_data.root idName=passingPresel outputFileName=EE
cmsRun fitterWithManyTemplatesPresel.py isMC=False inputFileName=TnPTree_data.root idName=passingPresel outputFileName=EB
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingPresel-EB.root -d PhotonToRECO -b -m passingPresel -p
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingPresel-EE.root -d PhotonToRECO -b -m passingPresel -p
```
- IDMVA
```
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/getTemplatesFromMC.py --input TnPTree_mc_newID.root -o templateIDMVAEB.root -d PhotonToRECO --idprobe passingIDMVAEB --var1Name probe_Pho_r9 --var1Bins 0,0.85,1 --var2Name probe_sc_abseta --var2Bins 0,1.479 --addProbeCut "passingPresel==1"
(write down MC efficiencies printed out by the script)
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/getTemplatesFromMC.py --input TnPTree_mc_newID.root -o templateIDMVAEE.root -d PhotonToRECO --idprobe passingIDMVAEE --var1Name probe_Pho_r9 --var1Bins 0,0.9,1 --var2Name probe_sc_abseta --var2Bins 1.479,2.5 --addProbeCut "passingPresel==1"
(write down MC efficiencies printed out by the script)
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitIDMVAEB.py --idLabel passingIDMVA -t templateIDMVAEB.root --var1Bins 0,0.85,1 --var2Bins 0,1.479  --passBkgPdf "RooExponential:backgroundPass(mass, aPass[-0.01,-1.,0.])" --failBkgPdf "RooExponential:backgroundFail(mass, aFail[-0.01,-1.,0.])"
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitIDMVAEE.py --idLabel passingIDMVA -t templateIDMVAEE.root --var1Bins 0,0.90,1 --var2Bins 1.479,2.5 --passBkgPdf "RooExponential:backgroundPass(mass, aPass[-0.01,-1.,0.])" --failBkgPdf "RooExponential:backgroundFail(mass, aFail[-0.01,-1.,0.])"
```
Before running the following commands check the binning defined in the fitterWithManyTemplates.py configuration (and include the appropriate include for EB or EE).
```
cmsRun fitterWithManyTemplates.py isMC=False inputFileName=TnP_data.root outputFileName=EB idName=passingIDMVA
cmsRun fitterWithManyTemplates.py isMC=False inputFileName=TnP_data.root outputFileName=EE idName=passingIDMVA
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingIDMVA-EB.root -d PhotonToRECO -b -m passingIDMVA -p 
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingIDMVA-EE.root -d PhotonToRECO -b -m passingIDMVA -p 
```

- Systematics 

Instructions to run systematics are shown only for preselection since there is no systematic error for the preselection cut on the IDMVA.

  - Running parametric fit for MC (Signal PDF systematic)
```
cmsRun fitter.py outputFileName=Systematics
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-mc-passingPresel-Systematics.root -d PhotonToRECO -b -m MCtruth_passingPresel
python getParamFromMCFit.py --input efficiency-mc-passingPresel-Systematics.root -d PhotonToRECO/MCtruth_passingPresel -p alphaF,alphaP,nF,nP --var1Name probe_Pho_r9 --var1Bins 0,0.85,1 --var2Name probe_sc_abseta --var2Bins 0,1.479
python getParamFromMCFit.py --input efficiency-mc-passingPresel-Systematics.root -d PhotonToRECO/MCtruth_passingPresel -p alphaF,alphaP,nF,nP --var1Name probe_Pho_r9 --var1Bins 0,0.9,1 --var2Name probe_sc_abseta --var2Bins 1.479,2.5
```
Modify the signal PDF in fitter.py fixing the parameters dervide from MC with the previous two commands and run the fit again.
```
cmsRun fitter.py isMC=False inputFileName=TnP_data.root outputFileName=SignalSyst
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingPresel-SignalSyst.root -d PhotonToRECO -b -m passingPresel
```
  - Changing background PDF 
```
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitPreselBkgSystEB.py --idLabel passingPresel -t templatePreselBkgSystEB.root --var1Bins 0,0.85,1 --var2Bins 0,1.479
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitPreselBkgSystEE.py --idLabel passingPresel -t templatePreselBkgSystEE.root --var1Bins 0,0.9,1 --var2Bins 1.479,2.5
```
Include the proper "commonFit..." file in fitterWithManyTemplatesPresel.py (and check again the binning):
```
cmsRun fitterWithManyTemplatesPresel.py isMC=False inputFileName=TnPTree_data.root idName=passingPresel outputFileName=BkgSystEB
cmsRun fitterWithManyTemplatesPresel.py isMC=False inputFileName=TnPTree_data.root idName=passingPresel outputFileName=BkgSystEE
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingPresel-BkgSystEB.root -d PhotonToRECO -b -m passingPresel -p
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingPresel-BkgSystEE.root -d PhotonToRECO -b -m passingPresel -p
```
Other systematics can be evaluated in similar fashion (tag cuts: estimated varying pT and id cuts on the tag, R9 reweighing: you can add a weight to the ntuple correspoinding to electron->photon R9 distribution and use this additional weight in the measurement i.e. using PhysicsTools/TagAndProbe/test/utilities/normalizeMCWeights.py script)
Beware that in some cases you may need to rerun the ntuplization step.

- OTHER CHECKS

You need to rerun the ntuplization step after having changed the photonIDMVA cut on the probes (e.g. 0.2, in ../python/treeMakerOptionsPhotons_cfi.py module process.goodPhotonProbesIDMVA)

  - IDMVA efficiency vs nVtx
```
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/getTemplatesFromMC.py --input TnP_mc_0p2.root -o templateID02Nvtx.root -d PhotonToRECO --idprobe passingIDMVA --var1Name event_nPV --var1Bins 0,5,9,13,17,100 --var2Name probe_sc_abseta --var2Bins 0,2.5 
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitIDMVA0p2Nvtx.py --idLabel passingIDMVA -t templateID02Nvtx.root --var1Bins 0,5,9,13,17,100 --var2Bins 0,2.5
cmsRun fitterWithManyTemplates.py isMC=False inputFileName=TnP_data.root idName=passingIDMVA outputFileName=IDMVA0p2_nvtx
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingIDMVA-IDMVA0p2_nvtx.root -d PhotonToRECO -b -m passingIDMVA
```
  - IDMVA efficiency vs pT
```
python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/getTemplatesFromMC.py --input TnP_mc_0p2.root -o templateID02pt.root -d PhotonToRECO --idprobe passingIDMVA --var1Name probe_Pho_et --var1Bins 20,30,40,50,60,1000 --var2Name probe_sc_abseta --var2Bins 0,2.5 python ../../../PhysicsTools/TagAndProbe/test/MCTemplates/makeConfigForTemplates.py -o ../python/commonFitIDMVA0p2pt.py --idLabel passingIDMVA -t templateID02pt.root --var1Bins 20,30,40,50,60,1000 --var2Bins 0,2.5
cmsRun fitterWithManyTemplates.py isMC=False inputFileName=TnP_data_0p2.root idName=passingIDMVA outputFileName=IDMVA0p2_et
python ../../../PhysicsTools/TagAndProbe/test/utilities/dumpPlotFromEffFile.py --input efficiency-data-passingIDMVA-IDMVA0p2_et.root -d PhotonToRECO -b -m passingIDMVA
```

#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.Systematics.SystematicDumperDefaultVariables import defaultVariables,minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables
import os

process = cms.Process("FLASHggTHQ")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' 
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )



from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process)
#modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

systlabels = [""]
phosystlabels = []
metsystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels = []

from flashgg.MetaData.JobConfig import customize
customize.options.register('doFiducial',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doFiducial'
                           )
customize.options.register('acceptance',
                           'NONE',
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.string,
                           'acceptance'
                           )
customize.options.register('doSystematics',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doSystematics'
                           )
customize.options.register('doPdfWeights',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doPdfWeights'
                           )
customize.options.register('dumpTrees',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'dumpTrees'
                           )
customize.options.register('dumpWorkspace',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'dumpWorkspace'
                           )


print "Printing defaults"
print 'doFiducial '+str(customize.doFiducial)
print 'acceptance '+str(customize.acceptance)
# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()
print "Printing options"
print 'doFiducial '+str(customize.doFiducial)
print 'acceptance '+str(customize.acceptance)

if customize.doFiducial:
    import flashgg.Systematics.fiducialCrossSectionsCustomize as fc
    fc.leadCut = 1./3.
    fc.subLeadCut = 1./4.
    fc.isoCut = 10.
    fc.etaCut = 2.5
    matchCut = "leadingPhoton.hasMatchedGenPhoton() && subLeadingPhoton.hasMatchedGenPhoton()"
    phoIDcut = '(leadingView().phoIdMvaWrtChosenVtx() >0.320 && subLeadingView().phoIdMvaWrtChosenVtx() >0.320)'
    accCut   = fc.getAccRecoCut()
    
    print process.flashggPreselectedDiPhotons.cut

    if customize.acceptance == 'IN':
        process.flashggPreselectedDiPhotons.cut = cms.string(str(process.flashggPreselectedDiPhotons.cut)[12:-2] +' && '+ str(matchCut)+ ' && '+ str(phoIDcut) +' && ' + str(accCut))

    if customize.acceptance == 'OUT':
        process.flashggPreselectedDiPhotons.cut = cms.string(str(process.flashggPreselectedDiPhotons.cut)[12:-2] +' && '+ str(matchCut)+ ' && '+ str(phoIDcut) +' && !' + str(accCut))
        
    if customize.acceptance == 'NONE':
        process.flashggPreselectedDiPhotons.cut = cms.string(str(process.flashggPreselectedDiPhotons.cut)[12:-2] +' && '+ str(phoIDcut))
    print "Here we print the preslection cut"
    print process.flashggPreselectedDiPhotons.cut

process.load("flashgg/Taggers/flashggTagSequence_cfi")
print 'here we print the tag sequence before'
print process.flashggTagSequence
if customize.doFiducial:
    from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet,massSearchReplaceAnyInputTag
    process.flashggTagSequence.remove(process.flashggVBFTag)
    process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
    #haven't tested VH tags with fiducial cross-section measurement yet
    process.flashggTagSequence.remove(process.flashggVHEtTag)
    process.flashggTagSequence.remove(process.flashggVHLooseTag)
    process.flashggTagSequence.remove(process.flashggVHTightTag)
    process.flashggTagSequence.remove(process.flashggVHMetTag)
    process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
    process.flashggTagSequence.remove(process.flashggVHHadronicTag)
    process.flashggTagSequence.replace(process.flashggUntagged, process.flashggSigmaMoMpToMTag)


print 'here we print the tag sequence after'
print process.flashggTagSequence

if customize.doFiducial:
    print 'we do fiducial and we change tagsorter'
    process.flashggTagSorter.TagPriorityRanges = cms.VPSet(     cms.PSet(TagName = cms.InputTag('flashggSigmaMoMpToMTag')) )



print "customize.processId:",customize.processId
# load appropriate scale and smearing bins here
# systematics customization scripts will take care of adjusting flashggDiPhotonSystematics
#process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")

# Or use the official tool instead
useEGMTools(process)

import flashgg.Taggers.THQLeptonicTagVariables as var
'''
variablesToUse = defaultVariables

variablesToUse.append("topMass :=getTopMass()")

variablesToUse.append("fwdJetEta :=getFwdJet().eta")

variablesToUse.append("LeptonPt :=getLeptonP4().pt")
variablesToUse.append("LeptonType :=getLeptonType()")
variablesToUse.append("LeptonCharge :=getLeptonCharge()")

variablesToUse.append("nJets :=jets().size()")
variablesToUse.append("nBJets :=bJets().size()")

variablesToUse.append("FoxWolf :=getFoxWolframMoment_ONE()")
variablesToUse.append("Aplanarity :=getAplanarity()")

variablesToUse.append("MVA :=thqleptonicMvaRes()")

variablesToUse.append("MET :=getMET()")
variablesToUse.append("METPhi :=getMET_Phi()")

variablesToUse.append("prompt_pho_1     := diPhoton.leadingPhoton.genMatchType()")
variablesToUse.append("prompt_pho_2     := diPhoton.subLeadingPhoton.genMatchType()")


variablesToUse.append("MuonsPt := getMuonsPt()" )
'''
variablesToUse = defaultVariables + var.vtx_variables + var.dipho_variables + var.photon_variables + var.lepton_variables + var.jet_variables + var.thqmva_variables
print "-------------------------------------------------"
print "--- Variables to be dumped, including systematic weights ---"
print variablesToUse
print "------------------------------------------------------------"


###### Dumper section

from flashgg.MetaData.samples_utils import SamplesManager

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/THQ_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_100016/0000/myMicroAODOutputFile_2.root"
                                                               )
                             )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"))

process.extraDumpers = cms.Sequence()
from flashgg.Taggers.tagsDumpers_cfi import *
#process.load("flashgg.Taggers.globalVariables_cff")

process.thqLeptonicTagDumper = createTagDumper("THQLeptonicTag")
process.thqLeptonicTagDumper.dumpTrees = customize.dumpTrees
process.thqLeptonicTagDumper.dumpWorkspace = customize.dumpWorkspace
process.thqLeptonicTagDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"
#process.thqLeptonicTagDumper.globalVariables = process.globalVariables

import flashgg.Taggers.dumperConfigTools as cfgTools
cfgTools.addCategories(process.thqLeptonicTagDumper,
                       ## categories definition
                       [("all","1",0)
                        ],                       
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       variables=variablesToUse,
                       histograms=[]
                       )

# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
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



if( not hasattr(process,"options") ): process.options = cms.untracked.PSet()
process.options.allowUnscheduled = cms.untracked.bool(True)


process.flashggTHQLeptonicTag.ElectronTag = "flashggElectronSystematics"
process.flashggTHQLeptonicTag.MuonTag = "flashggMuonSystematics"
process.p = cms.Path(process.dataRequirements*
                     process.genFilter*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     #process.flashggMets*
                     #process.flashggMetSystematics*
                     #process.flashggMuonSystematics*process.flashggElectronSystematics*
                     process.flashggUnpackedJets*
                     #process.jetSystematicsSequence*
                     process.flashggTagSequence*
                     #*process.systematicsTagSequences*
                     #process.flashggSystTagMerger*
                     #process.finalFilter*
                     process.thqLeptonicTagDumper)



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

# Detailed tag interpretation information printout (blinded)
process.flashggTagSorter.StoreOtherTagInfo = True
process.flashggTagSorter.BlindedSelectionPrintout = True

#### BELOW HERE IS MOSTLY DEBUGGING STUFF

#####################################################################
## Memory and timing, n.b. igprof is very complicated to interpret ##
##################################################################### 

#from Validation.Performance.TimeMemoryInfo import customise as TimeMemoryCustomize
#TimeMemoryCustomize(process)
#process.MessageLogger.cerr.threshold = 'WARNING'

#process.load("IgTools.IgProf.IgProfTrigger")
#process.igprof.reportEventInterval     = cms.untracked.int32(250)
#process.igprof.reportToFileAtBeginJob  = cms.untracked.string("|gzip -c>igprof.begin-job.gz")
#process.igprof.reportToFileAtEvent     = cms.untracked.string("|gzip -c>igprof.%I.%E.%L.%R.event.gz")
#process.p += process.igprof

################################
## Dump merged tags to screen ##
################################

#process.load("flashgg/Taggers/flashggTagTester_cfi")
#process.flashggTagTester.TagSorter = cms.InputTag("flashggSystTagMerger")
#process.flashggTagTester.ExpectMultiples = cms.untracked.bool(True)
#process.p += process.flashggTagTester

############################################
## Additional details on tag sorter steps ##
############################################

#process.flashggTagSorter.Debug = True

##############
## Dump EDM ##
##############

#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('CustomizeWillChangeThisAnyway.root'),
#                               outputCommands = cms.untracked.vstring('keep *') # dump everything! small tests only!
#                               )
#process.e = cms.EndPath(process.out)

############################
## Dump the output Python ##
############################
#print process.dumpPython()
#processDumpFile = open('processDump.py', 'w')
#print >> processDumpFile, process.dumpPython()

# set default options if needed
customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",1.00e+3)
# call the customization
customize(process)

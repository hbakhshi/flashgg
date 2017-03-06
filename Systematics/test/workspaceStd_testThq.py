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
process.load("flashgg.Taggers.flashggTagSequence_cfi")
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

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
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doSystematics'
                           )
customize.options.register('doPdfWeights',
                           False,
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
customize.options.register('doMuFilter',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doMuFilter'
                           )


from flashgg.MetaData.JobConfig import customize
customize.parse()



print process.flashggTagSequence




print "customize.processId:",customize.processId
# load appropriate scale and smearing bins here
# systematics customization scripts will take care of adjusting flashggDiPhotonSystematics
#process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")

# Or use the official tool instead
useEGMTools(process)


###### Dumper section

from flashgg.MetaData.samples_utils import SamplesManager

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
        "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/THQ_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_100016/0000/myMicroAODOutputFile_2.root"
        #"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016B-23Sep2016-v3/161114_162631/0000/myMicroAODOutputFile_122.root"
        #"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/THW_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_100138/0000/myMicroAODOutputFile_13.root"
        )
                             )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"))

process.extraDumpers = cms.Sequence()

import flashgg.Taggers.dumperConfigTools as cfgTools

import flashgg.Taggers.THQLeptonicTagVariables as var
variablesToUse = minimalVariables + defaultVariables + var.vtx_variables + var.dipho_variables + var.photon_variables + var.lepton_variables + var.jet_variables + var.thqmva_variables
print "-------------------------------------------------"
print "--- Variables to be dumped, including systematic weights ---"
print variablesToUse
print "------------------------------------------------------------"


from flashgg.Taggers.tagsDumpers_cfi import *

process.thqLeptonicTagDumper = createTagDumper("THQLeptonicTag")
process.thqLeptonicTagDumper.dumpTrees = customize.dumpTrees
process.thqLeptonicTagDumper.dumpWorkspace = customize.dumpWorkspace
process.thqLeptonicTagDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategories(process.thqLeptonicTagDumper,
                       ## categories definition
                       [("all","1",0)
                        ],                       
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       variables=variablesToUse,
                       histograms=[]
                       )


##### TAGDUMPER
process.load("flashgg.Taggers.diphotonTagDumper_cfi") ##  import diphotonTagDumper 
process.tagsDumper.className = "DiPhotonTagDumper"
process.tagsDumper.src = "flashggSystTagMerger"
#process.tagsDumper.src = "flashggTagSystematics"
process.tagsDumper.processId = "test"
process.tagsDumper.dumpTrees = customize.dumpTrees
process.tagsDumper.dumpWorkspace = customize.dumpWorkspace
process.tagsDumper.dumpHistos = False
process.tagsDumper.quietRooFit = True
process.tagsDumper.nameTemplate = cms.untracked.string("$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL")

tagList=[
    ["THQLeptonicTag",0]
    ]

definedSysts=set()
process.tagsDumper.classifierCfg.remap=cms.untracked.VPSet()
for tag in tagList: 
  tagName=tag[0]
  tagCats=tag[1]
  # remap return value of class-based classifier
  process.tagsDumper.classifierCfg.remap.append( cms.untracked.PSet( src=cms.untracked.string("flashgg%s"%tagName), dst=cms.untracked.string(tagName) ) )
  for systlabel in systlabels:
      if not systlabel in definedSysts:
          # the cut corresponding to the systematics can be defined just once
          cutstring = "hasSyst(\"%s\") "%(systlabel)
          definedSysts.add(systlabel)
      else:
          cutstring = None
      if systlabel == "":
          currentVariables = variablesToUse
      else:
          currentVariables = systematicVariables
      
      isBinnedOnly = (systlabel !=  "")
      if ( customize.doPdfWeights or customize.doSystematics ) and ( (customize.datasetName() and customize.datasetName().count("HToGG")) or customize.processId.count("h_") or customize.processId.count("vbf_") ) and (systlabel ==  ""):
          print "Signal MC central value, so dumping PDF weights"
          dumpPdfWeights = True
          nPdfWeights = 60
          nAlphaSWeights = 2
          nScaleWeights = 9
      else:
          print "Data, background MC, or non-central value, or no systematics: no PDF weights"
          dumpPdfWeights = False
          nPdfWeights = -1
          nAlphaSWeights = -1
          nScaleWeights = -1
      
      cfgTools.addCategory(process.tagsDumper,
                           systlabel,
                           classname=tagName,
                           cutbased=cutstring,
                           subcats=tagCats, 
                           variables=currentVariables,
                           histograms=minimalHistograms,
                           binnedOnly=isBinnedOnly,
                           dumpPdfWeights=dumpPdfWeights,
                           nPdfWeights=nPdfWeights,
                           nAlphaSWeights=nAlphaSWeights,
                           nScaleWeights=nScaleWeights
                           )


# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
# Bad Muon filter
process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
process.badGlobalMuonTagger.muons = cms.InputTag("flashggSelectedMuons")
process.cloneGlobalMuonTagger.muons = cms.InputTag("flashggSelectedMuons")

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter
        if customize.doMuFilter:
            process.dataRequirements += process.noBadGlobalMuons

process.genFilter = cms.Sequence()

if( not hasattr(process,"options") ): process.options = cms.untracked.PSet()
process.options.allowUnscheduled = cms.untracked.bool(True)


process.p = cms.Path(process.dataRequirements*
                     process.genFilter*
                     #process.flashggUpdatedIdMVADiPhotons*
                     #process.thqLeptonicTagDumper*
                     process.tagsDumper)



printSystematicInfo(process)

customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",1.00e+3)
# call the customization
customize(process)

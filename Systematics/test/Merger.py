import json
import os
import shutil
import sys
import subprocess

from ROOT import TFile, TTree, TH1, gDirectory, TList, TDirectory, gROOT, TObject, TString

class Merger:
    def __init__(self, name):
        self.Name = name 
        self.Jobs = []

    def AddJob(self, j):
        self.Jobs.append( j )

    def readKeys(self , directory):
        """Return a list of objects in directory that inherit from tree or histo. """

        if not directory.InheritsFrom("TDirectory"):
            return []
        selKeys = [key for key in directory.GetListOfKeys() if key.ReadObj().InheritsFrom("TH1") or key.ReadObj().InheritsFrom("TTree") or key.ReadObj().InheritsFrom("TDirectory")]
        ret = {}
        for k in selKeys:
            kcycle = k.GetCycle()
            kname = k.GetName()

            lastCycle = -1
            if kname in ret :
                lastCycle = ret[ kname ][0]
            if not (kcycle > lastCycle):
                continue
            elif (kcycle == lastCycle):
                print "%s has two similar cycle values %d and %d" % (kname , kcycle , lastCycle )

            ret[ kname ] = ( kcycle , k.ReadObj() )

        return [ ret[s][1] for s in ret ]


    def loop(self , directory):
        """Traverse directory recursively and return a list of (path, name) pairs of
        all objects inheriting from classname."""
        
        contents = []

        for d in self.readKeys(directory):
            if not d.InheritsFrom("TDirectory") :
                contents.append((directory.GetPath().split(':')[-1], d.GetName() ))
            else :
                contents += self.loop(d)

        return contents

    def hadd_workspaces(self):
        target = "WS_PreSel1Jet0B_" + self.Name + ".root"
        print ['hadd_workspaces', target ] + self.Jobs
        print subprocess.check_output(['hadd_workspaces', target ] + self.Jobs )

    def fhadd(self, force=False, verbose=False, slow=True):
        """ taken from https://root.cern.ch/phpBB3/viewtopic.php?t=14881
        This function will merge objects from a list of root files and write them    
        to a target root file. The target file is newly created and must not
        exist, or if -f ("force") is given, must not be one of the source files.
        
        IMPORTANT: It is required that all files have the same content!

        Fast but memory hungry alternative to ROOT's hadd.
        
        Arguments:

        target -- name of the target root file
        sources -- list of source root files
        classname -- restrict merging to objects inheriting from classname
        force -- overwrite target file if exists
        """

        target = self.Name + ".root"
        sources = self.Jobs

        TH1.AddDirectory(False)
        # check if target file exists and exit if it does and not in force mode
        if not force and os.path.exists(target):
            raise RuntimeError("target file %s exists" % target)

        # open the target file
        print "fhadd Target file:", target
        outfile = TFile(target, "RECREATE")

        # open the seed file - contents is looked up from here
        seedfilename = sources[0]
        print "fhadd Source file 1", seedfilename
        seedfile = TFile(seedfilename)

        # get contents of seed file
        print "looping over seed file"
        contents = self.loop(seedfile)
        print "done %d objects are ready to be merged" % len(contents)
        if( verbose ):
            for c in contents:
                print c
                

        # open remaining files
        otherfiles = []
        for n, f in enumerate(sources[1:]):
            #print "fhadd Source file %d: %s" % (n+2, f)
            otherfiles.append(TFile(f))

        
        cut = "" #(CMS_hgg_mass > 100 && CMS_hgg_mass < 180) && (diphoMVA > -0.4) && (n_loose_ele == 1 || n_LooseMu25 == 1) && (MET_pt > 30)"
        # loop over contents and merge objects from other files to seed file objects
        for n, (path, hname) in enumerate(contents):

            #print "fhadd Target object: %s" % os.path.join(path, hname)
            obj_path = os.path.join(path, hname)
            obj_ = seedfile.Get(obj_path[1:])

            outfile.cd('/')
            # create target directory structure
            for d in path.split('/')[1:]:
                directory = gDirectory.GetDirectory(d)
                if not directory:
                    gDirectory.mkdir(d).cd()
                else:
                    gDirectory.cd(d)
            obj = None
            IsTree = False

            #if "sigma" in obj_.GetName():
            #    continue

            if obj_.InheritsFrom("TTree"):
                #obj = obj_.CloneTree()
                obj = obj_.CopyTree( cut )
                IsTree = True
            else:
                continue
                obj = obj_.Clone()

            # merge objects
            l = TList()
            for o in [of.Get(obj_path[1:]) for of in otherfiles]:
                if IsTree :
                    l.Add( o.CopyTree( cut ) )
                else :
                    l.Add(o)
            obj.Merge(l)

            # delete objects if in slow mode
            if slow:
                #print "Deleting %d object(s)", l.GetEntries()
                l.Delete()

            # write object to target
            obj.Write(obj.GetName(), TObject.kOverwrite)

        print "Writing and closing file"

        # let ROOT forget about open files - prevents deletion of TKeys
        for f in [outfile, seedfile]+otherfiles:
            gROOT.GetListOfFiles().Remove(f);

        outfile.Write()
        outfile.Close()

        for f in [seedfile]+otherfiles:
            f.Close()

        print "fhadd completed"

jobs = None
# with open("bkg_jobs_3/task_config.json", "r" ) as cfin:
#     task_config = json.load(cfin)
#     jobs = task_config["jobs"]

with open("sig_jobs_6/task_config.json", "r" ) as cfin:
    task_config = json.load(cfin)
    jobs = task_config["jobs"]

# print len(jobs)
# task_config = None
#with open("data_jobs_forKirill_1/task_config.json", "r" ) as cfin:
#    task_config = json.load(cfin)
#with open("data_jobs_forKirill_2/task_config.json", "r" ) as cfin:
#    task_config = json.load(cfin)
#with open("data_jobs_forKirill_3/task_config.json", "r" ) as cfin:
#    task_config.update( json.load(cfin) )
#    jobs = task_config["jobs"]

# print len(jobs)

# FOR DATA
# dataset_names = [ "Run2016B-03Feb2017",
#                   "Run2016C-03Feb2017",
#                   "Run2016D-03Feb2017",
#                   "Run2016E-03Feb2017",
#                   "Run2016F-03Feb2017",
#                   "Run2016G-03Feb2017",
#                   "Run2016H-03Feb2017_ver2",
#                   "Run2016H-03Feb2017_ver3"
#                   ]

#FOR SIGNALS
dataset_names = [ "output_THQ_M125",
                  "output_THW_M125"
                 ]
#dataset_names = [ "TTH" , "THQ" , "THW" , "VH" , "GGH" , "VBF" ]
#dataset_names = [ #"GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf",
                  #"GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8" ,
                  #"QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8" ,
		  #"QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8",
		  #"QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8",
                  #"DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8" ,
                  #"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" ,
                  #"TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8",
                  #"TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8",
                  #"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
                  #"ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
                  #"TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8",
                  #"TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
                  #"TTGJets",
                  #"TTGG"
                  #]

datasets = { ds:Merger(ds) for ds in dataset_names }

for job in jobs:
    cmd, args, outfile_, nsub, ret, batchId = job

    outfile = str(outfile_)
    if not os.path.exists(outfile):
        # if "THQ" in outfile or "THW" in outfile :
        #     print "Skipping the job"
        #     continue
        # else:
        print outfile, "doesn't exist"
        #continue

    ds = ""
    for a in dataset_names:
        if a in outfile:
            ds = a

    if ds == "" :
        #print "no dataset found for %s" % outfile
        continue

    datasets[ ds ].AddJob( outfile )


for ds_n in datasets:
    ds = datasets[ds_n]

    if len(ds.Jobs) == 0 :
        continue

    print "start merging %s, #files %d" % ( ds.Name, len(ds.Jobs) )
    

    ds.fhadd(slow=True)
    #ds.hadd_workspaces()

    print "Done"


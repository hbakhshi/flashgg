import json
import os
import shutil
import sys
import subprocess

from ROOT import TFile, TTree, TH1, gDirectory, TList, TDirectory, gROOT, TObject

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
        target = "WS_" + self.Name + ".root"
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
            if obj_.InheritsFrom("TTree"):
                obj = obj_.CloneTree()
            else:
                obj = obj_.Clone()

            # merge objects
            l = TList()
            for o in [of.Get(obj_path[1:]) for of in otherfiles]:
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

jobs = None
with open("sig_jobs_3/task_config.json", "r" ) as cfin:
    task_config = json.load(cfin)
    jobs = task_config["jobs"]


# FOR DATA
# dataset_names = [ "Run2016B-23Sep2016",
#                   "Run2016C-23Sep2016",
#                   "Run2016D-23Sep2016",
#                   "Run2016E-23Sep2016",
#                   "Run2016F-23Sep2016",
#                   "Run2016G-23Sep2016",
#                   "Run2016H-PromptReco-v2",
#                   "Run2016H-PromptReco-v3"
#                   ]

#FOR SIGNALS
dataset_names = [ "VBFHToGG_M125_13TeV_powheg_pythia8",
                  "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8",
                  "ttHToGG_M125_13TeV_powheg_pythia8_v2",
                  "THQ_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1",
                  "THW_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1"
                  ]

datasets = { ds:Merger(ds) for ds in dataset_names }

for job in jobs:
    cmd, args, outfile_, nsub, ret, batchId = job

    outfile = str(outfile_)
    ds = ""
    for a in dataset_names:
        if a in outfile:
            ds = a

    if ds == "" :
        print "no dataset found for %s" % outfile
        continue

    datasets[ ds ].AddJob( outfile )


for ds_n in datasets:
    ds = datasets[ds_n]
    print "start merging %s, #files %d" % ( ds.Name, len(ds.Jobs) )
    

    ds.fhadd()
    ds.hadd_workspaces()

    print "Done"


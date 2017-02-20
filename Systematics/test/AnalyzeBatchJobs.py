import json
import os
import shutil

from ROOT import TFile, TTree

jobs = None
with open("data_jobs_3/task_config.json", "r" ) as cfin:
    task_config = json.load(cfin)
    jobs = task_config["jobs"]

status = {}
for job in jobs:
    cmd, args, outfile, nsub, ret, batchId = job

    # if os.path.isdir( "LSFJOB_{0:d}".format(batchId[1]) ):
    #     os.rename( "LSFJOB_{0:d}/STDOUT".format( batchId[1] ) , batchId[0] + ".log" )
    #     os.rename( "LSFJOB_{0:d}".format( batchId[1] ) , "LSFJOB_{0:d}_".format( batchId[1] ) )

    # continue

    cp = False
    OutFileIsThere=True
    if not os.path.isfile( str(outfile) ):
        OutFileIsThere = False
        cp = True
        status[ batchId[0] ] = "no file"
    else:
        f = TFile.Open( str(outfile) )
        if not f:
            status[ batchId[0] ] = "wrong file"
            cp = True

        else:
            d = f.GetDirectory("thqLeptonicTagDumper")
            if not d:
                cp = True
                status[ batchId[0] ] = "no dir"
            f.Close()
        
        # if cp :
            # os.remove(  str(outfile) )
    if cp :
        shutil.copy2( str(batchId[0])+".sh" , "./runit2/" )
        if OutFileIsThere:
            shutil.move(  str(outfile) , "./runit2/" )
    
for j in status :
    print j
    print status[j]

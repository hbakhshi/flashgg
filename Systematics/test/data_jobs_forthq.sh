# NB this command is specific to the configuration at IC and is not gaurenteed elsewhere
#LM=/afs/cern.ch/work/s/sethzenz/fromscratch107/CMSSW_8_0_8_patch1/src/flashgg/Systematics/test/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt #6.26/fb
#LM=/afs/cern.ch/work/s/sethzenz/fromscratch107/CMSSW_8_0_8_patch1/src/flashgg/MetaData/work/jsons/Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt #9.17/fb
queue="1nd"
#LM=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
LM=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
outdir="/afs/cern.ch/user/h/hbakhshi/work/tHq/CMSSW_8_0_25/src/flashgg/Systematics/test/"
version="1"
fggRunJobs.py --load data_jobs_forthq.json -d $outdir/data_jobs_${version} -x cmsRun workspaceStd_testThq.py maxEvents=-1 -n 500 -q ${queue} -D -P useAAA=1 --no-use-tarball lumiMask=${LM}
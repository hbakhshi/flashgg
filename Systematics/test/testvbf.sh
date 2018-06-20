#!/bin/bash
export XRD_NETWORKSTACK=IPv4
scp -p lxplus014:/tmp/x509up_u12330 .
export X509_USER_PROXY=$PWD/x509up_u12330
WD=$PWD
echo
echo
echo
cd /afs/cern.ch/work/h/hbakhshi/tHq/CMSSW_8_0_26_patch1
eval $(scram runtime -sh)
cd $WD
mkdir /afs/cern.ch/work/h/hbakhshi/tHq/CMSSW_8_0_26_patch1/src/flashgg/Systematics/test//sig_jobs_1
echo "ls $X509_USER_PROXY"
ls $X509_USER_PROXY
cmsRun /afs/cern.ch/work/h/hbakhshi/tHq/CMSSW_8_0_26_patch1/src/flashgg/Systematics/test/workspaceStd_testThq.py maxEvents=1000 useAAA=1 puTarget=2.39e+05,8.38e+05,2.31e+06,3.12e+06,4.48e+06,6e+06,7e+06,1.29e+07,3.53e+07,7.87e+07,1.77e+08,3.6e+08,6.03e+08,8.77e+08,1.17e+09,1.49e+09,1.76e+09,1.94e+09,2.05e+09,2.1e+09,2.13e+09,2.15e+09,2.13e+09,2.06e+09,1.96e+09,1.84e+09,1.7e+09,1.55e+09,1.4e+09,1.24e+09,1.09e+09,9.37e+08,7.92e+08,6.57e+08,5.34e+08,4.27e+08,3.35e+08,2.58e+08,1.94e+08,1.42e+08,1.01e+08,6.9e+07,4.55e+07,2.88e+07,1.75e+07,1.02e+07,5.64e+06,2.99e+06,1.51e+06,7.32e+05,3.4e+05,1.53e+05,6.74e+04,3.05e+04,1.52e+04,8.98e+03,6.5e+03,5.43e+03,4.89e+03,4.52e+03,4.21e+03,3.91e+03,3.61e+03,3.32e+03,3.03e+03,2.75e+03,2.47e+03,2.21e+03,1.97e+03,1.74e+03,1.52e+03,1.32e+03,1.14e+03,983,839 campaign=RunIISummer16-2_4_1-25ns_Moriond17 processIdMap=/afs/cern.ch/work/h/hbakhshi/tHq/CMSSW_8_0_26_patch1/src/flashgg/Systematics/test/sig_jobs_1/config.json dataset=/VBFHToGG_M-125_13TeV_powheg_pythia8 outputFile=/afs/cern.ch/work/h/hbakhshi/tHq/CMSSW_8_0_26_patch1/src/flashgg/Systematics/test//sig_jobs_1/output_VBFHToGG_M-125_13TeV_powheg_pythia8.root nJobs=19 jobId=2
retval=$?
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        cp -pv $file /afs/cern.ch/work/h/hbakhshi/tHq/CMSSW_8_0_26_patch1/src/flashgg/Systematics/test/sig_jobs_1
        if [[ $? != 0 ]]; then
            errors="$errors $file($?)"
        fi
    done
    if [[ -n "$errors" ]]; then
       echo "Errors while staging files"
       echo "$errors"
       exit -2
    fi
fi

exit $retval


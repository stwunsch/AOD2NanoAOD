#!/bin/bash

# Define path for job directories and processes
BASE_PATH=/ceph/wunsch/cms_opendata_2012_jobs
mkdir -p $BASE_PATH
PROCESSES=( \
    SMHiggsToZZTo4L \
    GluGluToHToTauTau \
    VBF_HToTauTau \
    ZZTo2e2mu \
    ZZTo4mu \
    ZZTo4e \
    TTbar \
    W1JetsToLNu \
    W2JetsToLNu \
    W3JetsToLNu \
    DYJetsToLL_M-50 \
    Run2012B_DoubleMuParked \
    Run2012C_DoubleMuParked\
    )

# Create JDL files and job directories
for PROCESS in ${PROCESSES[@]}
do
    ./create_jdl.py $PROCESS $BASE_PATH/$PROCESS
done

# Submit jobs
THIS_PWD=$PWD
for PROCESS in ${PROCESSES[@]}
do
    cd $BASE_PATH/$PROCESS
    condor_submit job.jdl
    cd $THIS_PWD
done

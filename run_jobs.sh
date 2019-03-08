#!/bin/bash

# Define path for job directories
BASE_PATH=/ceph/wunsch/cms_opendata_2012_jobs_2019-03-08
mkdir -p $BASE_PATH

# Set processes
PROCESSES=( \
    #SMHiggsToZZTo4L \
    #ZZTo2e2mu \
    #ZZTo4mu \
    #ZZTo4e \
    GluGluToHToTauTau \
    VBF_HToTauTau \
    TTbar \
    W1JetsToLNu \
    W2JetsToLNu \
    W3JetsToLNu \
    DYJetsToLL \
    DY2JetsToLL \
    DY3JetsToLL \
    DY4JetsToLL \
    Run2012B_DoubleMuParked \
    Run2012C_DoubleMuParked \
    Run2012B_SingleMuParked \
    Run2012C_SingleMuParked \
    Run2012B_DoubleElectron \
    Run2012C_DoubleElectron \
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

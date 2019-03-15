#!/bin/bash

# Exit on error
set -e

echo "### Begin of job"

ID=$1
echo "ID:" $ID

PROCESS=$2
echo "Process:" $PROCESS

FILE=$3
echo "File:" $FILE

OUTPUT_DIR=/ceph/wunsch/opendata_files_2019-03-15
echo "Output directory:" $OUTPUT_DIR

CMSSW_BASE=/portal/ekpbms2/home/wunsch/workspace/ttbar/CMSSW_5_3_32
echo "CMSSW base:" $CMSSW_BASE

if [[ ${FILE} == *"Run2012"* ]]; then
    CONFIG=${CMSSW_BASE}/src/workspace/AOD2NanoAOD/configs/data_cfg.py
else
    CONFIG=${CMSSW_BASE}/src/workspace/AOD2NanoAOD/configs/simulation_cfg.py
fi
echo "CMSSW config:" $CONFIG

echo "Hostname:" `hostname`

echo "How am I?" `id`

echo "Where am I?" `pwd`

echo "### Start working"

# Make output directory
mkdir -p ${OUTPUT_DIR}/${PROCESS}

# Setup CMSSW
THIS_DIR=$PWD
cd $CMSSW_BASE
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $THIS_DIR

# Copy config file
cp $CONFIG cfg.py

# Get lumi mask for data files
mkdir -p data
wget https://raw.githubusercontent.com/stwunsch/AOD2NanoAOD/master/data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
mv Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt data/

# Modify CMSSW config to run only a single file
sed -i -e "s,^files =,files = ['"${FILE}"'] #,g" cfg.py
sed -i -e 's,^files.extend,#files.extend,g' cfg.py

# Run CMSSW config
cmsRun cfg.py

# Copy output file
cp output.root ${OUTPUT_DIR}/${PROCESS}/${PROCESS}_${ID}.root

echo "### End of job"

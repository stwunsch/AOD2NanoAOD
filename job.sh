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

OUTPUT_DIR=/storage/c/wunsch/cms_opendata_2012_jobs
echo "Output directory:" $OUTPUT_DIR

CMSSW_BASE=/portal/ekpbms2/home/wunsch/workspace/CMSSW_5_3_32
echo "Output directory:" $OUTPUT_DIR

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
cp $CMSSW_BASE/src/workspace/AOD2NanoAOD/configs/${PROCESS}_cfg.py .

# Modify CMSSW config to run only a single file
sed -i -e "s,^files =,files = ['"${FILE}"'] #,g" ${PROCESS}_cfg.py
sed -i -e 's,^files.extend,#files.extend,g' ${PROCESS}_cfg.py

# Rename output file in CMSSW config
sed -i -e 's,'${PROCESS}.root','${PROCESS}_${ID}.root',g' ${PROCESS}_cfg.py

# Run CMSSW config
cmsRun ${PROCESS}_cfg.py

# Copy output file
cp ${PROCESS}.root ${OUTPUT_DIR}/${PROCESS}/${PROCESS}_${ID}.root

echo "### End of job"

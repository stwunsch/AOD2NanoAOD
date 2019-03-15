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

EOS_HOME=/eos/home-s/swunsch
echo "EOS home:" $EOS_HOME

OUTPUT_DIR=${EOS_HOME}/opendata_files/
echo "Output directory:" $OUTPUT_DIR

CMSSW_BASE=${EOS_HOME}/opendata_cmssw/CMSSW_5_3_32
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

# Trigger auto mount of EOS
ls -la $EOS_HOME

# Make output directory
mkdir -p ${OUTPUT_DIR}/${PROCESS}

# Setup CMSSW
THIS_DIR=$PWD
cd $CMSSW_BASE
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $THIS_DIR

# Copy config file
mkdir -p configs/
cp $CONFIG configs/cfg_${ID}.py

# Modify CMSSW config to run only a single file
sed -i -e "s,^files =,files = ['"${FILE}"'] #,g" cfg_${ID}.py
sed -i -e 's,^files.extend,#files.extend,g' cfg_${ID}.py

# Modify CMSSW config to read lumi mask from EOS
sed -i -e 's,data/Cert,'${CMSSW_BASE}'/src/data/workspace/AOD2NanoAOD/Cert,g' cfg_${ID}.py

# Modify config to write output directly to EOS
sed -i -e 's,output.root,'${OUTPUT_DIR}/${PROCESS}/${PROCESS}_${ID}'.root,g' cfg_${ID}.py

# Run CMSSW config
cmsRun cfg_${ID}.py

echo "### End of job"

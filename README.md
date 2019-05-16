# Convert AOD of CMS Open Data to NanoAOD

Tool to convert AOD to NanoAOD file format

## Setup CMSSW

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_5_3_32
```

## Build module

```bash
cd CMSSW_5_3_32/src
cmsenv
mkdir workspace
cd workspace
git clone <THIS REPOSITORY>
cd AOD2NanoAOD
scram b -j8
```

## Test configuration locally

```bash
cmsRun configs/simulation_cfg.py
cmsRun configs/data_cfg.py
```

## Create jobs on lxplus

```bash
./submit_jobs.sh /path/to/job/directory
```

## Merge job files

```bash
./merge_jobs.py /path/to/job/outputs
```

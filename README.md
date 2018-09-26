# Convert AOD of CMS Open Data to NanoAOD

TODO: Add top-level description

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
git clone <THIS REPOSITORY> aod2nanoaod
cd aod2nanoaod
scram b -j8
```

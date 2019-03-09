#!/bin/bash

BASE_PATH=/ceph/wunsch/cms_opendata_2012_nanoaod_2019-03-08/

for FOLDER in $(ls $BASE_PATH | grep -v .root)
do
    ./merge_jobs.py ${BASE_PATH}/${FOLDER} &
done

wait

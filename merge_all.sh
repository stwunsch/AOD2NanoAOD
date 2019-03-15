#!/bin/bash

BASE_PATH=/ceph/wunsch/opendata_files_2019-03-15/

for FOLDER in $(ls $BASE_PATH | grep -v .root)
do
    ./merge_jobs.py ${BASE_PATH}/${FOLDER} &
done

wait

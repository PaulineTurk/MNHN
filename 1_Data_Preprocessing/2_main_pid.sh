#!/bin/bash

# bash 2_main_pid.sh > 2.out 2>&1 &

listFASTA=/home/pauline/Bureau/MNHN_RESULT_FINAL/1_DATA/Pfam_FASTA/
echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'

start=`date +%s`

for entry in "$listFASTA"/*
    do
    nohup python3 2_main_pid.py $entry &
    done

wait
end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"

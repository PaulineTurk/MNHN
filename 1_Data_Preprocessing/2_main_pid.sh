#!/bin/bash

listFASTA=/home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_FASTA/

start=`date +%s`

for entry in "$listFASTA"/*
    do
    nohup python3 2_main_pid.py $entry &
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"

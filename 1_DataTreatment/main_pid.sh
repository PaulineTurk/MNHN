#!/bin/bash

listFASTA=/Users/pauline/Desktop/EXP_DIMANCHE/FASTA/


start=`date +%s`

for entry in "$listFASTA"/*
    do
    nohup python3 main_pid.py $entry &
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"

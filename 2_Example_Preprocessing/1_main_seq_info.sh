#!/bin/bash

# bash 1_main_seq_info.sh > 1_TRAIN.out 2>&1 &
# bash 1_main_seq_info.sh > 1_TEST.out 2>&1 &


echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'
echo
echo "PID:" $$
echo

# train_test="TRAIN"
train_test="TEST"
listFASTA=/home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_$train_test/

echo "listFASTA = /home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_$train_test/"


start=`date +%s`

for entry in "$listFASTA"/*
    do
    nohup python3 1_main_seq_info.py $entry $train_test &
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"


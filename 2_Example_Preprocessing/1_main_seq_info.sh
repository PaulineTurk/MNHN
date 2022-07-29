#!/bin/bash

listFASTA=/home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_train/

echo "listFASTA = /home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_train/"
echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'

start=`date +%s`

for entry in "$listFASTA"/*
    do
    nohup python3 1_main_seq_info.py $entry &
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"


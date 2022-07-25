#!/bin/bash


#listFASTA=/Users/pauline/Desktop/MNHN_RESULT/1_DATA/Pfam_split/Pfam_train/
listFASTA=/Users/pauline/Desktop/MNHN_RESULT_MINI/1_DATA/Pfam_split/Pfam_train/


start=`date +%s`

for entry in "$listFASTA"/*
    do

    nohup python3 main_ex_save.py $entry &
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"


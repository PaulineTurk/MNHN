#!/bin/bash


listFASTA=/Users/pauline/Desktop/MNHN_RESULT_MINI/2_5_EXAMPLES/EXAMPLES_6


start=`date +%s`

for entry in "$listFASTA"/*
    do

    nohup python3 main_ex_selection.py $entry &
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"

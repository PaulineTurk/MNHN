#!/bin/bash


listFASTA=/Users/pauline/Desktop/MNHN_RESULT_MINI/2_5_EXAMPLES/EXAMPLES_6


start=`date +%s`

COUNTER=0

for entry in "$listFASTA"/*
    do
    if $COUNTER<100;
    then
    let COUNTER++
    nohup python3 main_ex_selection.py $entry &

    else
    wait
    COUNTER=0
    fi
    done

end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"

#!/bin/bash

# bash 2_main_test_brier_naive_bayes_multi.sh > 2.out 2>&1 &

echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'

start=`date +%s`


general_output=OUTPUT
mkdir -p $general_output

name_folder=OUTPUT/TEST_MULTI_$$
mkdir $name_folder

Direction=(ol or dl dr)


for j in ${Direction[@]};
    do
    nohup python3 2_main_test_brier_naive_bayes_multi.py $j > $name_folder/"$j".txt 2>&1 &
    done;

wait
end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"
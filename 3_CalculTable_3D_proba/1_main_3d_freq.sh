#!/bin/bash

# bash 1_main_3d_freq.sh > 1.out 2>&1 &

echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'
echo
echo "PID:" $$
echo
start=`date +%s`

general_output=OUTPUT
mkdir -p $general_output


name_folder=$general_output/3D_FREQ_6_40_50_$$
mkdir $name_folder

listDirection=("ol" "or" "dl" "dr")

for direction in ${listDirection[@]};
    do
    nohup python3 1_main_3d_freq.py $direction > $name_folder/"$direction".out 2>&1 &
    done;

end=`date +%s`
runtime=$((end-start))
echo "3D FREQ"
echo "DONE IN: $runtime s"

#!/bin/bash

# bash main_3d_count.sh > main_3d_count.out 2>&1 &

echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'
echo
echo "PID:" $$
echo

start=`date +%s`

name_folder=OUTPUT_CUBE_6_40_50_$$
mkdir $name_folder

listDirection=("ol" "or" "dl" "dr")

for direction in ${listDirection[@]};
    do
    nohup python3 main_3d_count.py $direction > $name_folder/"$direction".out 2>&1 &
    done;

end=`date +%s`
runtime=$((end-start))
echo
echo "3D COUNT"
echo "DONE IN: $runtime s"

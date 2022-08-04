#!/bin/bash

# bash 2_main_selection_pc_3d.sh > 2.out 2>&1 &

general_output=OUTPUT
mkdir -p $general_output

name_folder=OUTPUT/SELECTION_PC_3D_UNI_$$
mkdir $name_folder

Direction=(ol or dl dr)


for j in ${Direction[@]};
    do
    echo "TIME START"
    date +'%d/%m/%Y %H:%M:%S'
    echo
    echo "PID:" $$
    echo
    nohup python3 2_main_selection_pc_3d.py $j > $name_folder/"$j".txt 2>&1 &
    done;

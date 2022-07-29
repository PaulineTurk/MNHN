#!/bin/bash

name_folder=OUTPUT_SELECTION_PC_2D_$$
mkdir $name_folder

myPseudoCounter=(0 0.00001 0.0001 0.001 0.01 0.1 1 10)


for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))
    name_file="$pid_inf"_"$pid_sup"
    mkdir $name_folder/$name_file

    for j in ${myPseudoCounter[@]};
        do
        nohup python3 selection_pc_2d.py $pid_inf $pid_sup $j > $name_folder/$name_file/"$j".txt 2>&1 &
        done;

    done;


#!/bin/bash

name_folder=OUTPUT_SELECTION_PC_2D
mkdir $name_folder

myArray=(0 0.00001 0.0001 0.001 0.01 0.1 1 10)

pid_inf=90
pid_sup=100


name_file="$pid_inf"_"$pid_sup"
mkdir $name_folder/$name_file


for j in ${myArray[@]};
    do
    nohup python3 selection_pc_2d.py $pid_inf $pid_sup $j > $name_folder/$name_file/"$j".txt 2>&1 &
    done;


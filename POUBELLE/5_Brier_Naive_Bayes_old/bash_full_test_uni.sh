#!/bin/bash

name_folder=OUTPUT_UNI
mkdir $name_folder

myArray=("origine" "destination")

for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))
    name_file="$pid_inf"_"$pid_sup"
    mkdir $name_folder/$name_file


    for j in ${myArray[@]};
        do
        nohup python3 -u 4_main_brier.py $pid_inf $pid_sup uni $j > $name_folder/$name_file/3d_count_"$i"_"$j"_$$.txt 2>&1 &
        done;

    done
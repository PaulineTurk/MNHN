#!/bin/bash

name_folder=OUTPUT
mkdir $name_folder
for ((i=0; i<=9; i++)); 
    do
    pid_inf=$(( 10*$i ))
    pid_sup=$(( 10*$i + 10 ))
    name_file="$pid_inf"_"$pid_sup"
    nohup python3 -u 2_main_table_2d.py $pid_inf $pid_sup > $name_folder/2d_$name_file.txt 2>&1 &
    done

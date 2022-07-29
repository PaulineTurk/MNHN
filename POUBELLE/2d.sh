#!/bin/bash

name_folder=OUTPUT_2D
mkdir $name_folder

for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))
    nohup python3 main_table_2d.py $pid_inf $pid_sup > $name_folder/"$pid_inf"_"$pid_sup".txt 2>&1 &
    done
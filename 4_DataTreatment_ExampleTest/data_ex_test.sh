#!/bin/bash

name_folder=OUTPUT_EX_TEST_$$
mkdir $name_folder


for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))

    nohup python3 -u main_data_exemple_test.py $pid_inf $pid_sup > $name_folder/"$pid_inf"_"$pid_sup".txt 2>&1 &

    done
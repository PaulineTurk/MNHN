#!/bin/bash

name_folder=OUTPUT
mkdir $name_folder


for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))

    nohup python3 -u 4_main_data_exemple_test.py $pid_inf $pid_sup > $name_folder/data_ex_test_"$pid_inf"_"$pid_sup".txt 2>&1 &

    done
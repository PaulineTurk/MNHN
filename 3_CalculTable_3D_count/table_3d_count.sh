#!/bin/bash

name_folder=OUTPUT_$$
mkdir $name_folder

myArray=("o" "d")

for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))
    name_file="$pid_inf"_"$pid_sup"
    mkdir $name_folder/$name_file

    for ((i=-10; i<=10; i++));
        do
        if [ $i -ne 0 ]; then
            for j in ${myArray[@]};
                do
                nohup python3 -u main_table_3d_count.py $pid_inf $pid_sup $i $j > $name_folder/$name_file/"$i"_"$j".txt 2>&1 &
                done;
        fi;
        done;
    done
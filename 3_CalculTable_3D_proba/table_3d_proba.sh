#!/bin/bash


name_folder=OUTPUT
mkdir $name_folder

myArray=("origine" "destination")

for ((x=0; x<=9; x++)); 
    do
    pid_inf=$(( 10*$x ))
    pid_sup=$(( 10*$x + 10 ))
    name_file="$pid_inf"_"$pid_sup"
    mkdir $name_folder/$name_file             # éventuellement mkdir -p

    for ((i=-10; i<=10; i++));
        do
        if [ $i -ne 0 ]; then
            for j in ${myArray[@]};
                do
                nohup python3 -u 3_main_table_3d_proba.py $pid_inf $pid_sup $i $j > $name_folder/$name_file/3d_proba_"$i"_"$j".txt 2>&1 &
                done;
        fi;
        done;
    done
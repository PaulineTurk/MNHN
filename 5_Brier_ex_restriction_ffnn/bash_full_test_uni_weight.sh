#!/bin/bash

name_folder=OUTPUT_UNI
mkdir $name_folder

myArray=("origine" "destination")
weight=(0.001 0.01 0.1 1 10)
method=("uni" "multi")
 
pid_inf=60
pid_sup=70
name_file="$pid_inf"_"$pid_sup"
mkdir $name_folder/$name_file

for m in ${method[@]};
    do
    for w in ${weight[@]};
        do
        for j in ${myArray[@]};
            do
            nohup python3 -u main_brier.py $pid_inf $pid_sup $m $j $w> $name_folder/$name_file/"$j"_"$w"_"$m"_$$.txt 2>&1 &
            done;
        done;
    done;
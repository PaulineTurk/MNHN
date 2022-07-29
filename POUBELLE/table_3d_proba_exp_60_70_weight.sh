#!/bin/bash


name_folder=exp_60_70_weight_$$
mkdir $name_folder

myArray=("o" "d")

list_weight=(0.001 0.01 0.1 1 10)

pid_inf=60
pid_sup=70

for w in ${list_weight[@]};
    do
    for ((i=-10; i<=10; i++));
        do
        if [ $i -ne 0 ]; then
            for j in ${myArray[@]};
                do
                nohup python3 -u pseudo_compte.py $pid_inf $pid_sup $i $j $w > $name_folder/"$i"_"$j"_"$w".txt 2>&1 &
                done;
        fi;
        done;
    done
#!/bin/bash


name_folder=exp_60_70_pc_$$
mkdir $name_folder

myArray=("o" "d")



pid_inf=60
pid_sup=70


for ((i=-10; i<=10; i++));
    do
    if [ $i -ne 0 ]; then
        for j in ${myArray[@]};
            do
            nohup python3 -u main_pseudo_counter_3d.py $pid_inf $pid_sup $i $j $w > $name_folder/"$i"_"$j".txt 2>&1 &
            done;
    fi;
    done

#!/bin/bash

name_folder=OUTPUT_SELECTION_PC_3D_test_$$
mkdir $name_folder


myPseudoCounter=(0 0.00001 0.0001 0.001 0.01 0.1 1.0 10.0 100.0 1000.0 10000.0 100000.0 1000000.0 10000000.0 100000000.0 1000000000.0 10000000000.0)

# myMethod=('uni' 'multi')
# myRef=('o' 'd')

#for ((x=0; x<=9; x++)); 
    #do
    #pid_inf=$(( 10*$x ))
    #pid_sup=$(( 10*$x + 10 ))
    pid_inf=60
    pid_sup=70
    # name_file="$pid_inf"_"$pid_sup"
    # mkdir $name_folder/$name_file

    for j in ${myPseudoCounter[@]};
        do
        nohup python3 selection_pc_3d.py $pid_inf $pid_sup "uni" "o" $j > $name_folder/$name_file/o_uni_"$j".txt 2>&1 &
        done

    #done;


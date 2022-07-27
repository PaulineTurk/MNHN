#!/bin/bash

name_folder=OUTPUT_CUBE_6_40_50_$$
mkdir $name_folder

listDirection=("ol" "or" "dl" "dr")

for direction in ${listDirection[@]};
    do
    nohup python3 main_cube_new.py $direction > $name_folder/"$direction".out 2>&1 &
    done;

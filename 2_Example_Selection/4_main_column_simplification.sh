#!/bin/bash

# bash 4_main_column_simplification.sh > 4_main_column_simplification.out 2>&1 &


echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'
echo
echo "PID:" $$
echo

# COUNT EX PER SEED
start=`date +%s`
path_folder_examples=/home/pauline/Bureau/MNHN_RESULT/2_EXAMPLES/EXAMPLES_6_40_50/

echo $path_folder_examples

count=0
batch_size=99



for entry in "$path_folder_examples"/*;
    do  
        if [ $count -lt $batch_size ];
        then 
            nohup python3 4_main_column_simplification.py $entry &
            let count=count+1
            
        else
            nohup python3 4_main_column_simplification.py $entry &
            let count=count+1
            wait
            count=0
        fi
    done


end=`date +%s`
runtime=$((end-start))
echo "COLUMN REMOVED"
echo "DONE IN: $runtime s"
echo

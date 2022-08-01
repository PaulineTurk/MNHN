#!/bin/bash

# bash 3_main_ex_per_seed.sh > 3_TRAIN.out 2>&1 &
# bash 3_main_ex_per_seed.sh > 3_TEST.out 2>&1 &


echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'
echo
echo "PID:" $$
echo

# train_test="TRAIN"
train_test="TEST"

# COUNT EX PER SEED
start=`date +%s`
global_path=/home/pauline/Bureau/MNHN_RESULT/2_EXAMPLES_$train_test
folder_name=EXAMPLES_6_40_50
new_csv_name=num_ex

listCSV="$global_path"/"$folder_name"/
echo $listCSV

COUNTER=0

echo "id_seed,num_ex" > "$global_path"/"$new_csv_name".csv

for entry in "$listCSV"/*;
    do
        let COUNTER=COUNTER+1
        num_line=$(cat $entry | wc -l)
        num_ex=$(( $num_line - 1))
        id_seed=$(basename $entry .csv)
        echo "$id_seed,$num_ex" >> "$global_path"/"$new_csv_name".csv
    done

echo
echo "n_files: $COUNTER"
end=`date +%s`
runtime=$((end-start))
echo "COUNT EX PER SEED"
echo "DONE IN: $runtime s"
echo


# FRACTION EX PER SEED
start=`date +%s`
path_original=/home/pauline/Bureau/MNHN_RESULT/2_EXAMPLES_$train_test/num_ex.csv
path_target=/home/pauline/Bureau/MNHN_RESULT/2_EXAMPLES_$train_test/frac_ex.csv

echo $path_original
echo $path_target

nohup python3 3_main_ex_per_seed.py $path_original $path_target
end=`date +%s`
runtime=$((end-start))
echo "FRACTION EX PER SEED"
echo "DONE IN: $runtime s"
echo

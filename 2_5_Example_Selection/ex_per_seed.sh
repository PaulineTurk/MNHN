#!/bin/bash

# bash ex_per_seed.sh > ex_per_seed.out 2>&1 &

start=`date +%s`

global_path=/home/pauline/Bureau/MNHN_RESULT/2_5_EXAMPLES
folder_name=EXAMPLES_6_40_50
new_csv_name=num_ex

listCSV="$global_path"/"$folder_name"/


echo "PID:" $$

echo "id_seed,num_ex" > "$global_path"/"$new_csv_name".csv

for entry in "$listCSV"/*;
    do
        num_line=$(cat $entry | wc -l)
        num_ex=$(( $num_line - 1))
        id_seed=$(basename $entry .csv)
        echo "$id_seed,$num_ex" >> "$global_path"/"$new_csv_name".csv
    done

echo
end=`date +%s`
runtime=$((end-start))
echo "DONE IN: $runtime s"


# echo "$host, `date`, checkout,$Time_checkout" >> log.csv
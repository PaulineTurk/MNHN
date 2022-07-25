#!/bin/bash

start=`date +%s`

listFASTA=/home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_train/

count=0
batch_size=100
echo "PID:" $$
n_file=0

for entry in "$listFASTA"/*;
    do
        if [ $count -lt $batch_size ];
        then 
            nohup python3 main_ex_save.py $entry &
            let count=count+1
            
        else
            wait
            echo
            date +"%H:%M"
            n_file=$(($n_file + $count))
            echo "n_file:" $n_file
            count=0

        fi
    done

echo
date +"%H:%M"
n_file_final=$(($n_file + $count))
echo "n_file_final:" $n_file_final

end=`date +%s`
runtime=$((end-start))

echo
echo "DONE IN: $runtime s"

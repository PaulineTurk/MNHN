#!/bin/bash

start=`date +%s`

listFASTA=/home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_train/

echo "listFASTA = /home/pauline/Bureau/MNHN_RESULT/1_DATA/Pfam_split/Pfam_train/"
echo "TIME START"
date +'%d/%m/%Y %H:%M:%S'

count=0
batch_size=99
echo "PID:" $$
n_file=0
echo


for entry in "$listFASTA"/*;
    do  
        if [ $count -lt $batch_size ];
        then 
            nohup python3 2_main_ex_save.py $entry &
            let count=count+1
            
        else
            nohup python3 2_main_ex_save.py $entry &
            let count=count+1
            wait
            n_file=$(($n_file + $count))
            count=0
        fi
    done


n_file_final=$(($n_file + $count))
echo "n_file_final:" $n_file_final

end=`date +%s`
runtime=$((end-start))

echo
echo "DONE IN: $runtime s"

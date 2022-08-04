"""
DATA PRE-PROCESSING:
nohup python3 4_main_data_split.py > 4.out 2>&1 &
"""

# IMPORTS
import sys
from pathlib import Path
import time
from datetime import datetime

import treatment.split as split
import utils.folder as folder


file = Path(__file__).resolve()
sys.path.append(file.parents[0])



# PARAMETERS
TRAIN_PERCENTAGE = 90

## main folders
DATA =  f"{file.parents[2]}/MNHN_RESULT/1_DATA"
NAME_CLUSTER_FOLDER = "Pfam_nonRedondant"
NAME_SPLIT_DATA_FOLDER = "Pfam_split"





# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

start = time.time()

print("\nDATA SPLIT TRAIN/TEST:")

path_folder_fasta_nonRedondant = f"{DATA}/{NAME_CLUSTER_FOLDER}"


print(f"TRAIN_PERCENTAGE: {TRAIN_PERCENTAGE}")
path_folder_data_split = f"{DATA}/{NAME_SPLIT_DATA_FOLDER}"
folder.creat_folder(path_folder_data_split)
split.data_split(path_folder_fasta_nonRedondant,
            path_folder_data_split, TRAIN_PERCENTAGE, "Pfam_TRAIN", "Pfam_TEST")

end = time.time()
print(f"\nDONE IN {'{:_}'.format(round(end - start, 4))} s")

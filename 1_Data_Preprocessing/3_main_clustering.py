"""
CLUSTERING:
nohup python3 3_main_clustering.py > 3.out 2>&1 &
"""

# IMPORTS
import sys
from pathlib import Path
import time
from datetime import datetime


import treatment.redundancy as redundancy
import treatment.description as description


file = Path(__file__).resolve()
sys.path.append(file.parents[0])



# PARAMETERS
CLUSTERING_PID = 99
ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

## main folders
DATA =  f"{file.parents[2]}/MNHN_RESULT_FINAL/1_DATA"
NAME_FASTA_FOLDER_UPPER = "Pfam_Upper"
NAME_PID_FOLDER = "PID"
NAME_CLUSTER_FOLDER = "Pfam_nonRedondant"






# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

start = time.time()

print(f"CLUSTERING_PID: {CLUSTERING_PID}")

path_folder_fasta_upper = f"{DATA}/{NAME_FASTA_FOLDER_UPPER}"
path_folder_pid = f"{DATA}/{NAME_PID_FOLDER}"

path_folder_fasta_nonRedondant = f"{DATA}/{NAME_CLUSTER_FOLDER}"
redundancy.multi_non_redundancy_correction(path_folder_fasta_upper,
                        path_folder_fasta_nonRedondant,
                        path_folder_pid,
                        ALPHABET, CLUSTERING_PID)

path_character_percentage = f"{DATA}/character_cluster.npy"
path_character_included_percentage = f"{DATA}/character_included_cluster.npy"
description.data_count(path_folder_fasta_nonRedondant, ALPHABET,
                       path_character_percentage,
                       path_character_included_percentage)

description.bar_plot_data_description(path_folder_fasta_nonRedondant,
                                    path_character_percentage, "ALL CHARACTERS")
description.bar_plot_data_description(path_folder_fasta_nonRedondant,
                    path_character_included_percentage , "STANDARD AMINO-ACIDS")


end = time.time()
print(f"\nDONE IN {'{:_}'.format(round(end - start, 4))} s")

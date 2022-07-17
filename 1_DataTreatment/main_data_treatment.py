"""
DATA PRE-PROCESSING:
python3 main_data_treatment.py > data_processing.txt 2>&1
"""

# IMPORTS
import sys
import os
from pathlib import Path
import time

import treatment.stockholm as stockholm
import treatment.capitalizer as capitalizer
import treatment.pid as pid
import treatment.redundancy as redundancy
import treatment.split as split
import treatment.description as description
import utils.folder as folder


file = Path(__file__).resolve()
sys.path.append(file.parents[0])



# PARAMETERS
FOLDER_SOURCE =  f"{file.parents[2]}/MNHN_RESULT"
DATA =  f"{file.parents[2]}/MNHN_RESULT/1_DATA"
NAME_MULTI_STOCKHOLM_FILE = "Pfam-A.seed"
NAME_MONO_STOCKHOLM_FOLDER = "Pfam_Stockholm"
NAME_FASTA_FOLDER = "Pfam_FASTA"
NAME_FASTA_FOLDER_UPPER = "Pfam_Upper"
NAME_PID_FOLDER = "PID"
NAME_CLUSTER_FOLDER = "Pfam_nonRedondant"
NAME_SPLIT_DATA_FOLDER = "Pfam_split"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
CLUSTERING_PID = 99
TRAIN_PERCENTAGE = 50






print("_______________________________________________________________________")
print("                         DATA PRE-PROCESSING                           ")
print("             OF PROTEIC MULTI-SEQUENCE-ALIGNMENTS (MSA)                ")
print("_______________________________________________________________________")

print("DATA PROCESSING START")
start = time.time()

# general folder managment
if not os.path.exists(FOLDER_SOURCE):
    os.makedirs(FOLDER_SOURCE)
if not os.path.exists(DATA):
    os.makedirs(DATA)

print("\nMULTI --> MONO STOCKHOLM:")
path_file_multi_stockholm = f"{file.parents[2]}/{NAME_MULTI_STOCKHOLM_FILE}"
path_folder_mono_stockholm = f"{DATA}/{NAME_MONO_STOCKHOLM_FOLDER}"
stockholm.stockholm_separator(path_file_multi_stockholm, path_folder_mono_stockholm)


print("\nCONVERSION FROM STOCKHOLM FORMAT TO FASTA FORMAT:")
path_folder_fasta = f"{DATA}/{NAME_FASTA_FOLDER}"
stockholm.multi_stockholm_to_fasta(path_folder_mono_stockholm, path_folder_fasta)

path_character_percentage = f"{DATA}/character_stockholm.npy"
path_character_included_percentage = f"{DATA}/character_included_stockholm.npy"
description.data_count(path_folder_fasta, ALPHABET,
                       path_character_percentage,
                       path_character_included_percentage)

description.bar_plot_data_description(path_folder_fasta, 
                                    path_character_percentage, "ALL CHARACTERS")
description.bar_plot_data_description(path_folder_fasta, 
                    path_character_included_percentage , "STANDARD AMINO-ACIDS")


print("\nUPPER CASE:")
path_folder_fasta_upper = f"{DATA}/{NAME_FASTA_FOLDER_UPPER}"
capitalizer.multi_capitalization(path_folder_fasta, path_folder_fasta_upper)

# path_character_percentage = f"{DATA}/character_upper.npy"
# path_character_included_percentage = f"{DATA}/character_included_upper.npy"
# description.data_count(path_folder_fasta_upper, ALPHABET,
#                        path_character_percentage,
#                        path_character_included_percentage)

# description.bar_plot_data_description(path_folder_fasta_upper, 
#                                     path_character_percentage, "ALL CHARACTERS")
# description.bar_plot_data_description(path_folder_fasta_upper, 
#                     path_character_included_percentage , "STANDARD AMINO-ACIDS")


print("\nPID")
path_folder_pid = f"{DATA}/{NAME_PID_FOLDER}"
pid.save_pid(path_folder_fasta_upper, path_folder_pid, ALPHABET)


print("\nCLUSTERING")
print(f"CLUSTERING_PID: {CLUSTERING_PID}")
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



print("\nDATA SPLIT TRAIN/TEST:")
print(f"TRAIN_PERCENTAGE: {TRAIN_PERCENTAGE}")
path_folder_data_split = f"{DATA}/{NAME_SPLIT_DATA_FOLDER}"
folder.creat_folder(path_folder_data_split)
split.data_split(path_folder_fasta_nonRedondant, 
            path_folder_data_split, TRAIN_PERCENTAGE, "Pfam_train", "Pfam_test")

end = time.time()
print(f"\nDATA PROCESSING DONE IN {round(end - start, 2)} s")

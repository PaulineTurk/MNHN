"""
DATA PRE-PROCESSING
"""

# IMPORTS
import sys
from pathlib import Path

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

DATA =  f"{file.parents[2]}/MNHN_RESULT/1_DATA_mini"
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
print("            OF PROTEIC MULTI-SEQUENCE-ALIGNEMENTS (MSA)                ")
print("_______________________________________________________________________")


print("\nfrom MULTI Stockholm MSA to MONO Stockholm MSA")
path_file_multi_stockholm = f"{DATA}/{NAME_MULTI_STOCKHOLM_FILE}"
path_folder_mono_stockholm = f"{DATA}/{NAME_MONO_STOCKHOLM_FOLDER}"
stockholm.stockholm_separator(path_file_multi_stockholm, path_folder_mono_stockholm)


print("\nfrom Stockholm format to FASTA format")
path_folder_fasta = f"{DATA}/{NAME_FASTA_FOLDER}"
stockholm.multi_stockholm_to_fasta(path_folder_mono_stockholm, path_folder_fasta)
# residu_count, total_residu, character_count, total_character = description.data_count(path_folder_fasta, ALPHABET)
# description.bar_plot_data_count(path_folder_fasta, residu_count, total_residu, "Standard amino-acid")
# description.bar_plot_data_count(path_folder_fasta, character_count, total_character , "Character")


print("\nUpper Case")
path_folder_fasta_upper = f"{DATA}/{NAME_FASTA_FOLDER_UPPER}"
capitalizer.multi_capitalization(path_folder_fasta, path_folder_fasta_upper)

residu_count, total_residu, character_count, total_character = description.data_count(path_folder_fasta_upper, ALPHABET)
description.bar_plot_data_count(path_folder_fasta_upper, residu_count, total_residu, "Standard amino-acid")
description.bar_plot_data_count(path_folder_fasta_upper, character_count, total_character , "Character")


print("\nPID")
path_folder_pid = f"{DATA}/{NAME_PID_FOLDER}"
pid.save_pid(path_folder_fasta_upper, path_folder_pid, ALPHABET)


print("\nClustering")
path_folder_fasta_nonRedondant = f"{DATA}/{NAME_CLUSTER_FOLDER}"
redundancy.multi_non_redundancy_correction(path_folder_fasta_upper, path_folder_fasta_nonRedondant, ALPHABET, CLUSTERING_PID)

residu_count, total_residu, character_count, total_character = description.data_count(path_folder_fasta_nonRedondant, ALPHABET)
description.bar_plot_data_count(path_folder_fasta_nonRedondant, residu_count, total_residu, "Standard amino-acid")
description.bar_plot_data_count(path_folder_fasta_nonRedondant, character_count, total_character , "Character")



print("\nData split: train/test")
path_folder_data_split = f"{DATA}/{NAME_SPLIT_DATA_FOLDER}"
folder.creat_folder(path_folder_data_split)
split.data_split(path_folder_fasta_nonRedondant, path_folder_data_split, TRAIN_PERCENTAGE, "Pfam_train", "Pfam_test")

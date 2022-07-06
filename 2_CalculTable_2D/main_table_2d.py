"""
NON-CONTEXTUAL INFORMATION
"""


# IMPORTS
import os
import argparse

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table2d.table2dfonction as table2dfonction
import utils.folder as folder



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
args = parser.parse_args()

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
DATA = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
NAME_FASTA_TRAIN_FOLDER = "Pfam_split/Pfam_train"
NAME_PID_FOLDER = "PID"
PSEUDO_COUNTER = 1


print("_______________________________________________________________________")
print("              ALPHABET CHARACTERS COUNT AND FREQUENCY                  ")
print("                      INDIVIDUALLY AND BY PAIR                         ")
print("_______________________________________________________________________")

# output folder management: global
IS_EXIST = os.path.exists(DATA_RESULT)
if not IS_EXIST:
    os.makedirs(DATA_RESULT)

path_folder_fasta = f"{DATA}/{NAME_FASTA_TRAIN_FOLDER}"
path_folder_pid = f"{DATA}/{NAME_PID_FOLDER}"


print(f"2D_TABLE {args.pid_inf},{args.pid_sup}")
# output folder management: pid range
path_folder_Result = f"{DATA_RESULT}/table_2d_{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER}"
path_folder_Result = folder.creat_folder(path_folder_Result)


# table_2d count
count_AA, nb_AA, count_coupleAA, nb_coupleAA = table2dfonction.multi_count_for_table_2d(
                                                path_folder_fasta, path_folder_pid,
                                                ALPHABET, args.pid_inf, args.pid_sup,
                                                path_folder_Result,
                                                PSEUDO_COUNTER)

# table_2d frequency
freq_AA, freq_coupleAA = table2dfonction.freq_for_table_2d(
                                         count_AA, nb_AA,
                                         count_coupleAA, nb_coupleAA,
                                         path_folder_Result)


# table_2d score (BLOSUM formula)
table_2d_score = table2dfonction.table_2d_score(freq_AA, freq_coupleAA,
                                                path_folder_Result, scale_factor = 2)


# table_2d probability
cond_proba = table2dfonction.table_2d_conditional_proba(freq_AA, freq_coupleAA, path_folder_Result)

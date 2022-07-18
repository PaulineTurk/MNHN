"""
NON-CONTEXTUAL INFORMATION:
python3 main_table_2d.py 60 70 > data_table_2d_$$.txt 2>&1
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
DATA = f"{file.parents[2]}/MNHN_RESULT_DRAFT/1_DATA"
DATA_RESULT_2D_COUNT = f"{file.parents[2]}/MNHN_RESULT_DRAFT/2_TABLE_2D"
NAME_FASTA_TRAIN_FOLDER = "Pfam_split/Pfam_train"
NAME_PID_FOLDER = "PID"
PSEUDO_COUNTER_INITIAL = 1
LIST_PSEUDO_COUNTER_2D = [0, 1, 10, 100, 1_000, 10_000, 100_000, 1_000_000]


print("_______________________________________________________________________")
print("              ALPHABET CHARACTERS COUNT AND FREQUENCY                  ")
print("                      INDIVIDUALLY AND BY PAIR                         ")
print("_______________________________________________________________________")

# output folder management: global

if not os.path.exists(DATA_RESULT_2D_COUNT):
    os.makedirs(DATA_RESULT_2D_COUNT)

path_folder_fasta = f"{DATA}/{NAME_FASTA_TRAIN_FOLDER}"
path_folder_pid = f"{DATA}/{NAME_PID_FOLDER}"


print(f"2D_TABLE {args.pid_inf},{args.pid_sup}")
print(f"PSEUDO_COUNTER_INITIAL: {PSEUDO_COUNTER_INITIAL}")

# output folder management: pid range
path_folder_Result = f"{DATA_RESULT_2D_COUNT}/{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER_INITIAL}"
path_folder_Result = folder.creat_folder(path_folder_Result)


print("")
print("_______________")
print(f"TABLE 2D COUNT")
print("_______________")
table2dfonction.multi_count_for_table_2d(path_folder_fasta, path_folder_pid,
                                         ALPHABET, args.pid_inf, args.pid_sup,
                                         path_folder_Result,
                                         PSEUDO_COUNTER_INITIAL)
path_count_AA = f"{path_folder_Result}/table_1d_count.npy"
path_count_AA_couple = f"{path_folder_Result}/table_2d_count.npy"



print("")
print("_______________")
print(f"TABLE 2D SCORE")
print("_______________")
table2dfonction.score(path_count_AA,
                      path_count_AA_couple,
                      path_folder_Result,
                      ALPHABET,
                      scale_factor=2)



for pseudo_counter_2d in LIST_PSEUDO_COUNTER_2D:
    print("")
    print("_________________________________")
    print(f"TABLE 2D PROBA")
    print(f"PSEUDO_COUNTER_2D: {pseudo_counter_2d}")
    print("_________________________________")
    table2dfonction.proba_conditional_weighted(path_count_AA, path_count_AA_couple,
                                               pseudo_counter_2d,
                                               ALPHABET,
                                               path_folder_Result)

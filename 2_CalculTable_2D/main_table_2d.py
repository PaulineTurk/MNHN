"""
NON-CONTEXTUAL INFORMATION:
python3 main_table_2d.py 90 100 > data_table_2d_90_100_$$.txt 2>&1
"""


# IMPORTS
import os
import argparse
import csv
import os.path

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
DATA_RESULT_2D_COUNT = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
NAME_FASTA_TRAIN_FOLDER = "Pfam_split/Pfam_train"
NAME_PID_FOLDER = "PID"
LIST_PSEUDO_COUNTER_2D = [0,
                          pow(10, -5),
                          pow(10, -4),
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          1,
                          10]


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

# output folder management: pid range
path_folder_Result = f"{DATA_RESULT_2D_COUNT}/{args.pid_inf}_{args.pid_sup}"
path_folder_Result = folder.creat_folder(path_folder_Result)


print("")
print("_______________")
print(f"TABLE 2D COUNT")
print("_______________")
percentage_not_evaluated = table2dfonction.multi_count_for_table_2d(
                                         path_folder_fasta,
                                         path_folder_pid,
                                         ALPHABET, args.pid_inf, args.pid_sup,
                                         path_folder_Result)

path_file_param_eval = f"{DATA_RESULT_2D_COUNT}/Parameters_Evaluated.csv"
file_exists = os.path.isfile(path_file_param_eval)

with open (path_file_param_eval, 'a', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header = ("pid_inf", "pid_sup", "percentage_not_evaluated")

    if not file_exists:
        writer.writerow(header)

    data = (args.pid_inf, args.pid_sup, percentage_not_evaluated)
    writer.writerow(data)


path_count_AA = f"{path_folder_Result}/count_1d.npy"
path_count_AA_couple = f"{path_folder_Result}/count_2d.npy"

print("")
print("_______________")
print(f"TABLE 2D FREQ")
print("_______________")

table2dfonction.freq(path_count_AA, path_count_AA_couple,
                     ALPHABET, path_folder_Result)


path_freq_AA = f"{path_folder_Result}/freq_1d.npy"
path_freq_AA_couple = f"{path_folder_Result}/freq_2d.npy"

print("")
print("_______________")
print(f"TABLE 2D SCORE")
print("_______________")
table2dfonction.score(path_freq_AA,
                      path_freq_AA_couple,
                      path_folder_Result,
                      ALPHABET,
                      scale_factor=2)



for pseudo_counter_2d in LIST_PSEUDO_COUNTER_2D:
    pseudo_counter_2d = float(pseudo_counter_2d)
    print("")
    print("_________________________________")
    print(f"TABLE 2D PROBA")
    print(f"PSEUDO_COUNTER_2D: {pseudo_counter_2d}")
    print("_________________________________")
    table2dfonction.proba_conditional_weighted(path_freq_AA, path_freq_AA_couple,
                                               pseudo_counter_2d,
                                               ALPHABET,
                                               path_folder_Result)

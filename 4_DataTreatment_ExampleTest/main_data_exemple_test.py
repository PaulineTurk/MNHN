"""
Preprocessing of data selection:
python3 main_data_exemple_test.py 90 100 > 90_100_$$.txt 2>&1
"""

# IMPORTS
import os
import os.path
import argparse
import csv

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import data_exemple_test.dataexampletestfonction as dataexampletestfonction
import utils.folder as folder

# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
args = parser.parse_args()


DATA = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
DATA_RESULT_EX_TEST = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST"
ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

NOM_FOLDER_FASTA_TEST = "Pfam_split/Pfam_test"
NOM_FOLDER_PID = "PID"



print("_______________________________________________________________________")
print("                           PRE-PROCESSING                              ")
print("                       EXAMPLE TEST SELECTION                          ")
print("_______________________________________________________________________")

print(f"PID_INF: {args.pid_inf}")
print(f"PID_SUP: {args.pid_sup}")


path_folder_seed = f"{DATA}/{NOM_FOLDER_FASTA_TEST}"   # seed test
path_folder_pid = f"{DATA}/{NOM_FOLDER_PID}"           # seed test pid

# FOLDER MANAGEMENT
path_review = f"{DATA_RESULT_EX_TEST}/REVIEW"
new_folder_dico = f"{DATA_RESULT_EX_TEST}/{args.pid_inf}_{args.pid_sup}"

list_path = [DATA_RESULT_EX_TEST, path_review, new_folder_dico]
for path in list_path:
    if not os.path.exists(path):
        os.makedirs(path)

# FOLDER CREATION
path_folder_dico_seq = folder.creat_folder(f"{new_folder_dico}/seq")
path_folder_dico_seed = folder.creat_folder(f"{new_folder_dico}/seed")

# DICO-SEED AND DICO-SEQ CALCULATION
n_valid_seed, n_non_info_seed, n_valid_aa_couple_global = dataexampletestfonction.multi_seeds_selection(
                                                            path_folder_seed, path_folder_pid,
                                                            args.pid_inf, args.pid_sup, ALPHABET,
                                                            path_folder_dico_seq, path_folder_dico_seed)

# REVIEW REGISTRATION

path_ex_test = f"{path_review}/example_test.csv"
file_exists = os.path.isfile(path_ex_test)

with open (path_ex_test, 'a', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header = ('pid_inf', 'pid_sup', 'n_valid_seed', 'n_non_info_seed', 'n_valid_aa_couple_global')

    if not file_exists:
        writer.writerow(header)

    data = args.pid_inf, args.pid_sup, n_valid_seed, n_non_info_seed, n_valid_aa_couple_global
    writer.writerow(data)

"""
Preprocessing of data selection
"""

# IMPORTS
import os
import argparse

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


path_folder_seed = f"{DATA}/{NOM_FOLDER_FASTA_TEST}"   # seed test
path_folder_pid = f"{DATA}/{NOM_FOLDER_PID}"           # seed test pid

# folder management: global
IS_EXIST = os.path.exists(DATA_RESULT_EX_TEST)
if not IS_EXIST:
    os.makedirs(DATA_RESULT_EX_TEST)


# folder management: per pid range
new_folder_dico = f"{DATA_RESULT_EX_TEST}/{args.pid_inf}_{args.pid_sup}"
IS_EXIST = os.path.exists(new_folder_dico)
if not IS_EXIST:
    os.makedirs(new_folder_dico)

# file creation
path_folder_dico_seq = folder.creat_folder(f"{new_folder_dico}/seq")
path_folder_dico_seed = folder.creat_folder(f"{new_folder_dico}/seed")

# dico_seed and dico_seq calculation
n_valid_seed, n_invalid_seed, n_valid_aa_couple = dataexampletestfonction.multi_seeds_selection(path_folder_seed, path_folder_pid,
                                                                                                args.pid_inf, args.pid_sup, ALPHABET,
                                                                                                path_folder_dico_seq, path_folder_dico_seed)

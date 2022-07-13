"""
NON-CONTEXTUAL INFORMATION WITH pc
"""


# IMPORTS
import os
import argparse
import numpy as np

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
parser.add_argument("pc", help="pc x freq(aa_origine) x freq(aa_destination)", type=float)
args = parser.parse_args()

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

DATA_RESULT_2D = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
PSEUDO_COUNTER = 1


print("_______________________________________________________________________")
print("                          2D pc WEIGHT                             ")
print("_______________________________________________________________________")

print(f"2D_TABLE {args.pid_inf},{args.pid_sup}")
print(f"PSEUDO-COMPTE: {args.pc}")

# output folder management: global
if not os.path.exists(DATA_RESULT_2D):
    os.makedirs(DATA_RESULT_2D)

# output folder management: pid range
path_folder_Result = f"{DATA_RESULT_2D}/{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER}"
if not os.path.exists(path_folder_Result):
    os.makedirs(path_folder_Result)

# table_2d count
# count_AA, nb_AA, count_coupleAA, nb_coupleAA = table2dfonction.multi_count_for_table_2d(
#                                                 path_folder_fasta, path_folder_pid,
#                                                 ALPHABET, args.pid_inf, args.pid_sup,
#                                                 path_folder_Result,
#                                                 PSEUDO_COUNTER)

# table_2d frequency_pc
path_count = f"{path_folder_Result}/count.npy"
path_freqAA = f"{path_folder_Result}/freq.npy"
path_Result = f"{path_folder_Result}/count_pc_{args.pc}.npy"
dico_count_2D_pc = table2dfonction.count_2d_pc(path_count, path_freqAA, args.pc, path_Result)
data = np.load(path_Result, allow_pickle=True).item()
print(table2dfonction.min_2D(data))

# # table_2d score (BLOSUM formula)
# table_2d_score = table2dfonction.table_2d_score(freq_AA, freq_coupleAA,
#                                                 path_folder_Result, scale_factor = 2)


# # table_2d probability
# cond_proba = table2dfonction.table_2d_conditional_proba(freq_AA, freq_coupleAA, path_folder_Result)

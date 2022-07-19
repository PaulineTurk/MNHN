"""
CONTEXTUAL INFORMATION
"""

# IMPORTS
import os
import argparse

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table3d_count.table3dcountfonction as table3dcountfonction



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("neighbour_relative_position",
                    help="position du voisin par rapport la position de référence",
                    type=int)
parser.add_argument("ref_seq", help="'o' for origine OR 'd' for destination", type=str)
args = parser.parse_args()


DATA = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
DATA_RESULT_3D = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

NAME_FASTA_TRAIN_FOLDER = "Pfam_split/Pfam_train"
NAME_PID_FOLDER = "PID"



print("_______________________________________________________________________")
print("                           3D_TABLE COUNT                              ")
print("                          BY RANGE OF PID                              ")
print("_______________________________________________________________________")

path_folder_fasta      = f"{DATA}/{NAME_FASTA_TRAIN_FOLDER}"
path_folder_pid        = f"{DATA}/{NAME_PID_FOLDER}"

# folder management: global 3D
if not os.path.exists(DATA_RESULT_3D):
    os.makedirs(DATA_RESULT_3D)

# folder management: global 3D_COUNT
DATA_RESULT_3D_COUNT = f"{DATA_RESULT_3D}/COUNT"
if not os.path.exists(DATA_RESULT_3D_COUNT):
    os.makedirs(DATA_RESULT_3D_COUNT)

# folder management: per PID range
path_table_3d_count = f"{DATA_RESULT_3D_COUNT}/{args.pid_inf}_{args.pid_sup}"
if not os.path.exists(path_table_3d_count):
    os.makedirs(path_table_3d_count)

# table_3d_count
print(f"PID_INF: {args.pid_inf}")
print(f"PID_SUP: {args.pid_sup}")
print(f"CONTEXT RELATIVE POSITION: {args.neighbour_relative_position}")
print(f"REFERENCE: {args.ref_seq}")
table3dcountfonction.multi_table_3d_count(path_folder_fasta, path_folder_pid,
                                        args.neighbour_relative_position, args.ref_seq,
                                        ALPHABET, args.pid_inf, args.pid_sup,
                                        path_table_3d_count)

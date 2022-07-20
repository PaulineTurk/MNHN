# IMPORTS
import sys
import argparse
import time
from pathlib import Path
import os
import numpy as np

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table3d_proba.table3dprobafonctioncorrection as table3dprobafonction



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("relative_c",
                    help="position du voisin par rapport la position de référence",
                    type=int)
parser.add_argument("ref_seq", help="'o' OR 'd'", type=str)

args = parser.parse_args()

file = Path(__file__).resolve()
sys.path.append(file.parents[1])


ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
DATA_2D = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
DATA_3D_COUNT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/COUNT"
DATA_3D_FREQ = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/FREQ"
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA"
PSEUDO_COUNTER_2D = pow(10,-2)
LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -3),
                          pow(10, -1),
                          pow(10, 0),
                          pow(10, 1),
                          pow(10, 5)]


table_1d = np.load(f"{DATA_2D}/{args.pid_inf}_{args.pid_sup}/freq_1d.npy", allow_pickle='TRUE').item()



for pseudo_counter_3d in LIST_PSEUDO_COUNTER_3D:
    pseudo_counter_3d = float(pseudo_counter_3d)
    print(f"PSEUDO_COUNTER_3D: {pseudo_counter_3d}")
    sum_cases = 0
    path_table_3d_proba = f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}/{pseudo_counter_3d}/{args.relative_c}_{args.ref_seq}.npy"
    table_proba_3d = np.load(path_table_3d_proba, allow_pickle='TRUE').item()

    for aa_c in ALPHABET:
        for aa_1 in ALPHABET:
            for aa_2 in ALPHABET:
                diff = abs(table_proba_3d[aa_1][aa_2][aa_c] - table_1d[aa_c])
                #print(f"{aa_1},{aa_2}: {diff}")
                sum_cases += diff
    print(sum_cases)
                


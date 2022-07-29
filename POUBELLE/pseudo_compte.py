"""
3D_TABLES 2D_WEIGHTED:
python3 pseudo_compte.py 60 70 1 o > 90_100_$$.txt 2>&1
"""

# IMPORTS
import sys
import argparse
from pathlib import Path
import os

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table3d_proba.table3dprobafonction as table3dprobafonction



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("relative_context_position",
                    help="position du voisin par rapport la position de référence",
                    type=int)
parser.add_argument("ref_seq", help="'o' OR 'd'", type=str)
parser.add_argument("weight", help="percentage of 2D weighted to 3D", type=float)
args = parser.parse_args()

file = Path(__file__).resolve()
sys.path.append(file.parents[1])


ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
DATA_3D_COUNT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/COUNT"
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA"
PSEUDO_COUNTER_2D = 1
PSEUDO_COUNTER_3D = 0


print("_______________________________________________________________________")
print("                           3D_TABLE PROBA                              ")
print("                     WEIGHTING and NORMALISATION                       ")
print("_______________________________________________________________________")


print("")
print("_________________")
print("    PARAMETERS   ")
print("_________________")

print(f"PID_INF: {args.pid_inf}")
print(f"PID_SUP: {args.pid_sup}")
print(f"RELATIVE_CONTEXT_POSITION: {args.relative_context_position}")
print(f"REFERENCE: {args.ref_seq}")
print(f"PSEUDO_COUNTER_2D: {PSEUDO_COUNTER_2D}")
print(f"PSEUDO_COUNTER_3D: {PSEUDO_COUNTER_3D}")
print(f"WEIGHT: {args.weight}")

# folder managment: global
list_data_folder_path = [DATA_3D_COUNT, DATA_3D_PROBA, DATA_2D_PROBA]
for path in list_data_folder_path:
    IS_EXIST = os.path.exists(path)
    if not IS_EXIST :
        os.makedirs(path)

# folder managment: per pid range
path_folder_save = f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER_3D}_{args.weight}"
IS_EXIST = os.path.exists(path_folder_save)
if not IS_EXIST :
    os.makedirs(path_folder_save)

print("")
print("_________________")
print("   CALCULATION   ")
print("_________________")
# table_3d_proba calculation, weighting and normalisation
table_3d_proba_w_n = table3dprobafonction.proba_normal_weighting(
                                            DATA_3D_COUNT,
                                            DATA_2D_PROBA,
                                            args.pid_inf, args.pid_sup,
                                            PSEUDO_COUNTER_2D, PSEUDO_COUNTER_3D,
                                            args.relative_context_position, args.ref_seq,
                                            args.weight, ALPHABET,
                                            DATA_3D_PROBA)
print("")
print("_________________________")
print("   CHECK NORMALISATION   ")
print("_________________________")
table3dprobafonction.sum_plate(table_3d_proba_w_n)

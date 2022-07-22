"""
3D_TABLES CONDITIONAL PROBABILITY WITH PSEUDO-COUNTERS:
python3 main_pseudo_counter_3d.py 90 100 1 o > 90_100_$$.txt 2>&1
"""

# IMPORTS
import sys
import argparse
import time
from pathlib import Path
import os

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table3d_proba.table3dprobafonctioncorrection as table3dprobafonction



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("relative_context_position",
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
                          pow(10, -5),
                          pow(10, -4),
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10, 0),
                          pow(10, 1),
                          pow(10, 2),
                          pow(10, 3),
                          pow(10, 4),
                          pow(10, 5),
                          pow(10, 6),
                          pow(10, 7),
                          pow(10, 8),
                          pow(10, 9),
                          pow(10, 10)]


print("_______________________________________________________________________")
print("                       3D_TABLE CONDITIONAL PROBA                      ")
print("                          WITH PSEUDO-COUNTERS                         ")
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


# folder managment: global
list_data_folder_path = [DATA_3D_COUNT,
                         DATA_3D_FREQ,
                         f"{DATA_3D_FREQ}/{args.pid_inf}_{args.pid_sup}",
                         DATA_3D_PROBA,
                         f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}",
                         DATA_2D]
for path in list_data_folder_path:
    if not os.path.exists(path):
        os.makedirs(path)

print("")
print("_________________")
print("  3D FREQUENCES  ")
print("_________________")
# table_3d_proba calculation, weighting and normalisation
path_table_3d_count =f"{DATA_3D_COUNT}/{args.pid_inf}_{args.pid_sup}/{args.relative_context_position}_{args.ref_seq}.npy"
path_save_table_3d_freq =f"{DATA_3D_FREQ}/{args.pid_inf}_{args.pid_sup}/{args.relative_context_position}_{args.ref_seq}.npy"

start_time = time.time()
table_3d_proba_w_n = table3dprobafonction.freq_3d(path_table_3d_count,
                                                      path_save_table_3d_freq,
                                                      ALPHABET)
end_time = time.time()
print(f"time: {round(end_time - start_time, 2)} s")


table3dprobafonction.sum_freq_3D(path_save_table_3d_freq, ALPHABET)



print("")
print("________________________________")
print("  3D PROBA WITH PSEUDO_COUNTER  ")
print("________________________________")

path_table_2d_proba_pc = f"{DATA_2D}/{args.pid_inf}_{args.pid_sup}/proba_{PSEUDO_COUNTER_2D}.npy"
path_table_1d_freq = f"{DATA_2D}/{args.pid_inf}_{args.pid_sup}/freq_1d.npy"

print(f"\nPATH_TABLE_3D_FREQ: {path_save_table_3d_freq}")
print(f"PATH_TABLE_2D_PROBA_PC: {path_table_2d_proba_pc}")
print(f"PATH_TABLE_1D_FREQ: {path_table_1d_freq}")



for pseudo_counter_3d in LIST_PSEUDO_COUNTER_3D:
    pseudo_counter_3d = float(pseudo_counter_3d)


    print(f"\nPSEUDO_COUNTER_3D: {pseudo_counter_3d}")
    
    # folder managment: per pid range and pseudo_counter_3d tested
    path_folder_table_3d_proba_pc = f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}/{pseudo_counter_3d}"
    if not os.path.exists(path_folder_table_3d_proba_pc):
        os.makedirs(path_folder_table_3d_proba_pc)

    path_file_table_3d_proba_pc = f"{path_folder_table_3d_proba_pc}/{args.relative_context_position}_{args.ref_seq}"

    start_time = time.time()
    table3dprobafonction.proba_3d(path_save_table_3d_freq,
                                  path_table_2d_proba_pc,
                                  path_table_1d_freq,
                                  pseudo_counter_3d,
                                  ALPHABET,
                                  path_file_table_3d_proba_pc)
    end_time = time.time()
    print(f"time: {round(end_time - start_time, 2)} s")
   
print("\nDONE")

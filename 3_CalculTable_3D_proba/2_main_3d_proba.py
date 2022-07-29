"""
3D_PROBA WEIGHTED
"""

# IMPORTS
import sys
import argparse
from pathlib import Path
import os
import time

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table3d_proba.function_3dproba as function_3dproba



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("direction",
                    help="'ol', 'or', 'dl', 'dr'", type=str)
args = parser.parse_args()

file = Path(__file__).resolve()
sys.path.append(file.parents[1])

L = 6
PID_INF = 40
PID_SUP = 50
ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

DATA_2D = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"
DATA_3D_COUNT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/COUNT/{L}_{PID_INF}_{PID_SUP}"
DATA_3D_PROBA_GLOBAL = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA"
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA/{L}_{PID_INF}_{PID_SUP}"
DATA_3D_FREQ = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/FREQ/{L}_{PID_INF}_{PID_SUP}"

PSEUDO_COUNTER_2D = pow(10, -2)

LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10, 0),
                          pow(10, 1),
                          pow(10, 2),
                          pow(10, 3)]


# PROGRAM

# folder managment
list_data_folder_path = [DATA_3D_PROBA, DATA_3D_PROBA_GLOBAL]
for path in list_data_folder_path:
    os.makedirs(path, exist_ok=True)



path_table_2d_proba_pc = f"{DATA_2D}/proba_{PSEUDO_COUNTER_2D}.npy"
path_table_1d_freq = f"{DATA_2D}/freq_1d.npy"
print(f"PATH_TABLE_2D_PROBA_PC: {path_table_2d_proba_pc}")
print(f"PATH_TABLE_1D_FREQ: {path_table_1d_freq}")

start_time = time.time()

for i in range(1, L+1):
    path_save_table_3d_freq = f"{DATA_3D_FREQ}/{args.direction}_{i}.npy"
    print(f"\nPATH_TABLE_3D_FREQ: {path_save_table_3d_freq}")


    for pseudo_counter_3d in LIST_PSEUDO_COUNTER_3D:
        pseudo_counter_3d = float(pseudo_counter_3d)

        print(f"\nPSEUDO_COUNTER_3D: {pseudo_counter_3d}")
        
        # folder managment
        path_folder_table_3d_proba_pc = f"{DATA_3D_PROBA}/{pseudo_counter_3d}"
        os.makedirs(path_folder_table_3d_proba_pc, exist_ok=True)

        path_file_table_3d_proba_pc = f"{path_folder_table_3d_proba_pc}/{args.direction}_{i}"

        function_3dproba.proba_3d(path_save_table_3d_freq,
                                  path_table_2d_proba_pc,
                                  path_table_1d_freq,
                                  pseudo_counter_3d,
                                  ALPHABET,
                                  path_file_table_3d_proba_pc)


end_time = time.time()
print(f"\nDONE IN: {round(end_time - start_time, 2)} s ")

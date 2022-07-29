"""
3D_FREQ
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

DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"
DATA_3D_COUNT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/COUNT/{L}_{PID_INF}_{PID_SUP}"
DATA_3D_FREQ_GLOBAL = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/FREQ"
DATA_3D_FREQ = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/FREQ/{L}_{PID_INF}_{PID_SUP}"



# PROGRAM

# folder managment
for path in [DATA_3D_FREQ, DATA_3D_FREQ_GLOBAL]:
    os.makedirs(path, exist_ok=True)

start_time = time.time()

for i in range(1, L+1):
    path_table_3d_count =f"{DATA_3D_COUNT}/{args.direction}_{i}.npy"
    path_save_table_3d_freq =f"{DATA_3D_FREQ}/{args.direction}_{i}.npy"


    function_3dproba.freq_3d(path_table_3d_count,
                            path_save_table_3d_freq,
                            ALPHABET)
    function_3dproba.sum_freq_3D(path_save_table_3d_freq, ALPHABET)

end_time = time.time()
print(f"\nDONE IN: {round(end_time - start_time, 2)} s")

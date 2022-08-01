"""
2D COUNT
nohup python3 1_main_2d_count.py > 1.out 2>&1 &
"""

# IMPORTS

import os
import os.path
import time
import csv
import numpy as np
from datetime import datetime

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])




# PARAMETERS
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)


PID_INF, PID_SUP = 40, 50
L = 6

DATA = f"{file.parents[2]}/MNHN_RESULT/2_EXAMPLES/EXAMPLES_{L}_{PID_INF}_{PID_SUP}"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


DATA_RESULT_GLOBAL = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D"
os.makedirs(DATA_RESULT_GLOBAL, exist_ok=True)
DATA_RESULT = f"{DATA_RESULT_GLOBAL}/{L}_{PID_INF}_{PID_SUP}"
os.makedirs(DATA_RESULT, exist_ok=True)


# PROGRAM

files = [x for x in Path(DATA).iterdir()]
n_files = len(files)


# initialisation of the count_2D
dict_2D = {}
for aa_origin in ALPHABET:
    dict_2D[aa_origin] = {}
    for aa_destination in ALPHABET:
        dict_2D[aa_origin][aa_destination] = 0


# count of the count_2D
counter = 0

start = time.time()
for file in files:
    with open(file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        print(file)
        for row in reader:
            dict_2D[row['aa_origin']][row['aa_destination']] += 1
    counter += 1
    print(round(100*counter/n_files, 2))

# save of the count_2D
np.save(f"{DATA_RESULT}/count_2d", dict_2D)


end = time.time()

print(f"DONE IN: {round(end - start, 2)} s")

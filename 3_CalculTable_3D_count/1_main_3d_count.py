"""
3D COUNT
"""

# IMPORTS

import os
import os.path
import argparse
import time
import csv
import numpy as np
from datetime import datetime

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table3d_count.table3dcountfonction as table3dcountfonction



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("direction",
                     help="'ol', 'or', 'dl', 'dr'", type=str)
args = parser.parse_args()

PID_INF, PID_SUP = 40, 50
L = 6

DATA = f"{file.parents[2]}/MNHN_RESULT/2_PRE_EXAMPLE_TRAIN/EXAMPLES_{L}_{PID_INF}_{PID_SUP}"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


DATA_RESULT_GLOBAL = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D"
DATA_RESULT_GLOBAL_COUNT = f"{DATA_RESULT_GLOBAL}/COUNT"
DATA_RESULT = f"{DATA_RESULT_GLOBAL}/COUNT/{L}_{PID_INF}_{PID_SUP}"

for path in [DATA_RESULT_GLOBAL, DATA_RESULT_GLOBAL_COUNT, DATA_RESULT]:
    os.makedirs(path, exist_ok=True)


# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

files = [x for x in Path(DATA).iterdir()]
n_files = len(files)

LIST_ABS_POSITION = [i for i in range(1, L+1)]


# initialisation of the count_3D
dict_dict_3D = {}
for position in LIST_ABS_POSITION:
    dict_dict_3D[(args.direction, position)] = table3dcountfonction.table_3d_count_initialisation(ALPHABET)


# count of the count_3D
counter = 0

start = time.time()
for file in files:
    with open(file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        print(file)
        for row in reader:
            for position in LIST_ABS_POSITION:
                dict_dict_3D[(args.direction, position)][row['aa_origin']][row['aa_destination']][row[f'aa_{args.direction}_{position}']] += 1
    counter += 1
    print(round(100*counter/n_files, 2))

# save of the count_3D
for position in LIST_ABS_POSITION:
    np.save(f"{DATA_RESULT}/{args.direction}_{position}", dict_dict_3D[(args.direction, position)])


end = time.time()

print(f"DONE IN: {round(end - start, 2)} s")

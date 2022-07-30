"""
SELECTION PSEUDO_COUNTER 3D - MULTI
"""

# IMPORTS
import os
import os.path
import numpy as np
import csv
import time
import argparse

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])


import brierNeighbour.brier as brier



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("direction", help="'ol', 'or', 'dl', 'dr'", type=str)
args = parser.parse_args()

L = 6
PID_INF = 40
PID_SUP = 50
# NAME_FILE_EXAMPLES = "EX_BRIER_TRAIN"
NAME_FILE_EXAMPLES = "EX_BRIER_TRAIN_1M"
NAME_EXPERIMENT = "EXP_1M_MULTI"

PSEUDO_COUNTER_2D = pow(10, -2)

LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10,  0),
                          pow(10,  1),
                          pow(10,  2),
                          pow(10,  3)]


# PATH FOR 2D_PROBA
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR 3D_PROBA
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLES = f"{file.parents[2]}/MNHN_RESULT/4_EXAMPLES/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/5_PC_3D_SELECTION"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]




# PROGRAM

# FOLDER MANAGMENT
os.makedirs(DATA_RESULT, exist_ok=True)
path_experiment = f"{DATA_RESULT}/{NAME_EXPERIMENT}"
os.makedirs(path_experiment, exist_ok=True)

## EXAMPLES
path_examples = f"{DATA_EXEMPLES}/{NAME_FILE_EXAMPLES}.csv"



start = time.time()

# LOAD: PROBA_2D_PC & 3D
# & INITIALISATION BRIER SCORE
dict_2d = np.load(f"{DATA_2D_PROBA}/proba_{PSEUDO_COUNTER_2D}.npy", allow_pickle='True').item()
dict_3d = {}
dict_score = {}

for position in range(0, L+1):
    dict_3d[position] = {}
    dict_score[position] = {}
    for pseudo_counter_3D in LIST_PSEUDO_COUNTER_3D:
        pseudo_counter_3D = float(pseudo_counter_3D)

        dict_score[position][pseudo_counter_3D] = []
        if position != 0:
            path_3d = f"{DATA_3D_PROBA}/{pseudo_counter_3D}/{args.direction}_{position}.npy"
            dict_3d[position][pseudo_counter_3D] = np.load(path_3d, allow_pickle='True').item()

            

# print(dict_score)

with open(path_examples, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        list_aa_context = []
        for position in range(0, L+1):
            aa_origin = row['aa_origin']
            aa_destination = row['aa_destination']
            if position != 0:
                list_aa_context.append(row[f'aa_{args.direction}_{position}'])

            for pseudo_counter_3D in LIST_PSEUDO_COUNTER_3D:
                pseudo_counter_3D = float(pseudo_counter_3D)

                # GET THE VECTOR OF PROBABILITY
                vect = brier.vecteur_from_table_3d_proba_multi(aa_origin,
                                                         list_aa_context,
                                                                 dict_2d,
                                                                 dict_3d,
                                                                position,
                                                       pseudo_counter_3D,
                                                                ALPHABET)

                brier_unit = brier.unit_brier_naive_bayes(vect, aa_destination, ALPHABET)
                dict_score[position][pseudo_counter_3D].append(brier_unit)

# print(dict_score)

np.save(f"{path_experiment}/SCORE_{args.direction}_multi", dict_score)


end = time.time()
print("")
print(f"DONE in {round(end-start, 4)} s")

"""
NON-CONTEXTUAL BRIER SCORE EXPERIMENT
"""

# IMPORTS
import os
import os.path
import numpy as np
import csv
import time

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])


import brierNeighbour.brier as brier



# PARAMETERS


L = 6
PID_INF = 40
PID_SUP = 50
#NAME_FILE_EXAMPLES = "EX_BRIER_TRAIN"
NAME_FILE_EXAMPLES = "EX_BRIER_TRAIN_1M"
NAME_EXPERIMENT = "EXP_1M"

LIST_PSEUDO_COUNTER_2D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10,  0),
                          pow(10,  1),
                          pow(10,  2),
                          pow(10,  3)]


# PATH FOR 2D_PROBA
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLES = f"{file.parents[2]}/MNHN_RESULT/4_EXAMPLES/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/5_PC_2D_SELECTION"

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

# LOAD: PROBA_2D_PC
# & INITIALISATION BRIER SCORE
dict_2d = {}
dict_score = {}

for pseudo_counter_2D in LIST_PSEUDO_COUNTER_2D:
    pseudo_counter_2D = float(pseudo_counter_2D)

    path_2d = f"{DATA_2D_PROBA}/proba_{pseudo_counter_2D}.npy"
    dict_2d[pseudo_counter_2D] = np.load(path_2d, allow_pickle='True').item()

    dict_score[pseudo_counter_2D] = []



with open(path_examples, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        for pseudo_counter_2D in LIST_PSEUDO_COUNTER_2D:
            pseudo_counter_2D = float(pseudo_counter_2D)

            # GET THE VECTOR OF PROBABILITY
            vect = []
            for aa in ALPHABET:
                vect.append(dict_2d[pseudo_counter_2D][row['aa_origin']][aa])
            brier_unit = brier.unit_brier_naive_bayes(vect, row['aa_destination'], ALPHABET)
            dict_score[pseudo_counter_2D].append(brier_unit)

# print(dict_score)

np.save(f"{path_experiment}/SCORE", dict_score)


end = time.time()
print("")
print(f"DONE in {round(end-start, 4)} s")

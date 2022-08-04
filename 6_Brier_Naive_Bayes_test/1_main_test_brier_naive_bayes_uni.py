"""
TEST BRIER NAIVE BAYES - UNI
"""

# IMPORTS
import os
import os.path
import numpy as np
import csv
import time
import argparse
from datetime import datetime

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
# TRAIN_TEST = "TRAIN
TRAIN_TEST = "TEST"


#NAME_FILE_EXAMPLES = "EX_BRIER_{TRAIN_TEST}"
NAME_FILE_EXAMPLES = f"EX_BRIER_{TRAIN_TEST}_1M"
NAME_EXPERIMENT = "EXP_1M_UNI"

PSEUDO_COUNTER_2D = pow(10, -2)
PSEUDO_COUNTER_3D = pow(10, -2)


# PATH FOR 2D_PROBA
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR 3D_PROBA
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLES = f"{file.parents[2]}/MNHN_RESULT/4_EXAMPLE_{TRAIN_TEST}/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/6_TEST_BRIER_NAIVE_BAYES"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]




# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)



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
    dict_score[position] = []
    if position != 0:
        path_3d = f"{DATA_3D_PROBA}/{PSEUDO_COUNTER_3D}/{args.direction}_{position}.npy"
        dict_3d[position][PSEUDO_COUNTER_3D] = np.load(path_3d, allow_pickle='True').item()

            

# print(dict_score)

with open(path_examples, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        for position in range(0, L+1):
            aa_origin = row['aa_origin']
            aa_destination = row['aa_destination']
            if position != 0:
                aa_context = row[f'aa_{args.direction}_{position}']
            else:
                aa_context = "None"


            # GET THE VECTOR OF PROBABILITY
            vect = brier.vecteur_from_table_3d_proba_uni(aa_origin,
                                                            aa_context,
                                                               dict_2d,
                                                               dict_3d,
                                                              position,
                                                     PSEUDO_COUNTER_3D,
                                                              ALPHABET)

            brier_unit = brier.unit_brier_naive_bayes(vect, aa_destination, ALPHABET)
            dict_score[position].append(brier_unit)

# print(dict_score)

np.save(f"{path_experiment}/SCORE_{args.direction}_uni", dict_score)


end = time.time()
print("")
print(f"DONE in {round(end-start, 4)} s")

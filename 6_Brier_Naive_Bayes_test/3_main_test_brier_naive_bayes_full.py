"""
TEST BRIER NAIVE BAYES - FULL
nohup python3 3_main_test_brier_naive_bayes_full.py > 3_10M.out 2>&1 &
"""

# IMPORTS
import os
import os.path
import numpy as np
import csv
import time
from datetime import datetime

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])


import brierNeighbour.brier as brier



# PARAMETERS
L = 6
PID_INF = 40
PID_SUP = 50


NAME_FILE_EXAMPLES = f"EX_BRIER_TEST_10M"
NAME_EXPERIMENT = "EXP_10M_FULL"


PSEUDO_COUNTER_2D = pow(10, -2)
LIST_PSEUDO_COUNTER_3D = [pow(10, 0)]



# PATH FOR 2D_PROBA
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR 3D_PROBA
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLES = f"{file.parents[2]}/MNHN_RESULT/4_EXAMPLE_TEST_10M/{L}_{PID_INF}_{PID_SUP}"
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/6_TEST_BRIER_NAIVE_BAYES"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]




# PROGRAM

# FOLDER MANAGMENT
os.makedirs(DATA_RESULT, exist_ok=True)
path_experiment = f"{DATA_RESULT}/{NAME_EXPERIMENT}"
os.makedirs(path_experiment, exist_ok=True)

## EXAMPLES
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

path_examples = f"{DATA_EXEMPLES}/{NAME_FILE_EXAMPLES}.csv"



start = time.time()

# LOAD: PROBA_2D_PC & 3D
# & INITIALISATION BRIER SCORE
dict_2d = np.load(f"{DATA_2D_PROBA}/proba_{PSEUDO_COUNTER_2D}.npy", allow_pickle='True').item()
dict_3d = {}
dict_score = {}

for direction in ["ol", "or", "dl", "dr"]:
    for position in range(0, L+1):
        dict_3d[(direction, position)] = {}
        dict_score[position] = {}
        for pseudo_counter_3D in LIST_PSEUDO_COUNTER_3D:
            pseudo_counter_3D = float(pseudo_counter_3D)

            dict_score[position][pseudo_counter_3D] = []
            if position != 0:
                path_3d = f"{DATA_3D_PROBA}/{pseudo_counter_3D}/{direction}_{position}.npy"
                dict_3d[(direction, position)][pseudo_counter_3D] = np.load(path_3d, allow_pickle='True').item()

            

# print(dict_score)

for pseudo_counter_3D in LIST_PSEUDO_COUNTER_3D:
    pseudo_counter_3D = float(pseudo_counter_3D)   # TEST PSEUDO-COUNTER-3D ALWAYS ON FIXED PSEUDO-COUNTER 2D and 3D

with open(path_examples, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        aa_origin = row['aa_origin']
        aa_destination = row['aa_destination']
        list_aa_context = []
        
        for position in range(0, L+1):
            if position != 0:
                for direction in ["ol", "or", "dl", "dr"]:
                    list_aa_context.append(row[f'aa_{direction}_{position}'])

            
            # GET THE VECTOR OF PROBABILITY
            vect = brier.vecteur_from_table_3d_proba_full(aa_origin,
                                                        list_aa_context,
                                                        dict_2d,
                                                        dict_3d,
                                                        position,
                                                        pseudo_counter_3D,
                                                        ALPHABET)

            brier_unit = brier.unit_brier_naive_bayes(vect, aa_destination, ALPHABET)
            dict_score[position][pseudo_counter_3D].append(brier_unit)

# print(dict_score)

np.save(f"{path_experiment}/SCORE_FULL", dict_score)


end = time.time()
print("")
print(f"DONE in {round(end-start, 4)} s")

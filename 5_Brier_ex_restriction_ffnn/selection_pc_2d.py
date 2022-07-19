"""
Brier Score Experiment:
nohup python3 main_brier_NEW.py 90 100 uni origine 0.001> test_brier_$$.txt 2>&1
"""

# IMPORTS
import os
import os.path
import argparse
import csv
import time

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import brierNeighbour.selection_example as selection_example
import brierNeighbour.brier as brier



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("pseudo_counter_2d",
                    help="pseudo_counter in 2D tables", type=float)
args = parser.parse_args()

# PATH FOR THE MSA ACCESS
DATA_MNHN = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
# PATH FOR 2D_PROBA
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLE_TEST = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST"
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/5_PC_2D_SELECTION"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
NAME_FASTA_TEST_FOLDER = "Pfam_split/Pfam_test"
N_TEST_PER_CONTEXT = 5
N_EXAMPLES_PER_TEST = 1_000_000


print("_______________________________________________________________________")
print("                     PSEUDO_COUNTER_2D SELECTION                       ")
print("_______________________________________________________________________")

  
print("_________________________________")
print("  Preparation of the experiment  ")
print("_________________________________")

list_context = [(0,0,0,0)]

print(f"PID_INF: {args.pid_inf}")
print(f"PID_SUP: {args.pid_sup}")
print(f"LIST CONTEXT TESTED: {list_context}")
print(f"N_TEST_PER_CONTEXT: {N_TEST_PER_CONTEXT}")
print(f"N_EXAMPLES_PER_TEST: {N_EXAMPLES_PER_TEST}")
print(f"PSEUDO_COUNTER_2D: {args.pseudo_counter_2d}")

list_nb_example = [N_EXAMPLES_PER_TEST]*N_TEST_PER_CONTEXT

# FOLDER MANAGMENT: global
range_pid_experiment = f"{DATA_RESULT}/Experiment_Brier_PC_2D_SELECTION"
list_path_folder = [DATA_RESULT, range_pid_experiment]
for path_folder in list_path_folder:
    if not os.path.exists(path_folder):
        os.makedirs(path_folder,  exist_ok=True)

# PATH TO LOAD:

## MSA
path_folder_seed = f"{DATA_MNHN}/{NAME_FASTA_TEST_FOLDER}"

## 2D_PROBA WITH PSEUDO_COUNTER_2D
path_table_2d_proba = f"{DATA_2D_PROBA}/{args.pid_inf}_{args.pid_sup}/proba_{args.pseudo_counter_2d}.npy"

## DATA PRE-PROECESSING DICO
path_folder_dico_seq = f"{DATA_EXEMPLE_TEST}/{args.pid_inf}_{args.pid_sup}/seq"
path_file_dico_seed_normal = f"{DATA_EXEMPLE_TEST}/{args.pid_inf}_{args.pid_sup}/seed/seed_normal.npy"


print("_________________________________")
print("            Experiment           ")
print("_________________________________")

start = time.time()

for context in list_context:
    print("")
    print("__________________________")
    print(f"CONTEXT: {context}")
    print("__________________________")

    context_ol, context_or, context_dl, context_dr = context

    test_number = 0
    total_list_unit_score_brier = []
    for nb_exemple_test in list_nb_example:
        test_number += 1
        print("")
        print(f"TEST {test_number}/{N_TEST_PER_CONTEXT}")
        dico_example = selection_example.example_number_per_seed(path_file_dico_seed_normal, N_EXAMPLES_PER_TEST)
        list_example = selection_example.multi_random_example_selection(path_folder_seed, dico_example,
                                                                        path_folder_dico_seq,
                                                                        context_ol, context_or, context_dl, context_dr,
                                                                        ALPHABET)
    
        score_brier, list_unit_score_brier, nb_example = brier.brier_score_no_context(
                                                                    list_example,
                                                                    ALPHABET,
                                                                    path_table_2d_proba)


        # REVIEW REGISTRATION

        path_selection_pc_2d = f"{range_pid_experiment}/{args.pid_inf}_{args.pid_sup}.csv"
        with open (path_selection_pc_2d, 'a', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            header = ('pseudo_counter_2d', 'score_brier', 'nb_example')   # prob with the header

            if not os.path.isfile(path_selection_pc_2d):
                writer.writerow(header)

            data = args.pseudo_counter_2d, score_brier, nb_example
            writer.writerow(data)


end = time.time()
print("")
print(f"DONE in {round(end-start, 4)} s")

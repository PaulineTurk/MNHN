"""
Brier Score Experiment:
bash bash_selection_pc_3d.sh
"""

# IMPORTS
import os
import argparse
import csv
import time
import numpy as np

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
parser.add_argument("method", help="'uni' OR 'multi'", type=str)
parser.add_argument("reference", help="'o' for origine OR 'd' for destination", type=str)
parser.add_argument("pseudo_counter_3d", help="pseudo_counter_3d", type=float)
args = parser.parse_args()

# PATH FOR THE MSA ACCESS
DATA_MNHN = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLE_TEST = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST"
# PATH FOR 2D_PROBA
DATA_2D_PROBA = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
# PATH FOR 3D_PROBA
DATA_3D_PROBA = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D/PROBA"
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/5_PC_3D_SELECTION"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
NAME_FASTA_TEST_FOLDER = "Pfam_split/Pfam_test"
MAX_RELATIVE_POSITION = 5
N_TEST_PER_CONTEXT = 1
N_EXAMPLES_PER_TEST = 100_000
PSEUDO_COUNTER_2D = pow(10, -2)


print("_______________________________________________________________________")
print("                     PSEUDO_COUNTER_3D SELECTION                       ")
print("_______________________________________________________________________")

  
print("_________________________________")
print("  Preparation of the experiment  ")
print("_________________________________")

if args.reference == "o":
    list_ol          = [(i, 0, 0, 0) for i in range(1,MAX_RELATIVE_POSITION + 1)[::-1]]
    list_non_context = [(0, 0, 0, 0)]
    list_or          = [(0, i, 0, 0) for i in range(1,MAX_RELATIVE_POSITION + 1)]

    list_context = list_ol  + list_non_context + list_or


if args.reference == "d":
    list_dl          = [(0, 0, i, 0) for i in range(1,MAX_RELATIVE_POSITION + 1)[::-1]]
    list_non_context = [(0, 0, 0, 0)]
    list_dr          = [(0, 0, 0, i) for i in range(1,MAX_RELATIVE_POSITION + 1)]

    list_context = list_dl  + list_non_context + list_dr

print(f"PID_INF: {args.pid_inf}")
print(f"PID_SUP: {args.pid_sup}")
print(f"METHOD: {args.method}")
print(f"LIST CONTEXT TESTED: {list_context}")
print(f"MAX_RELATIVE_POSITION: {MAX_RELATIVE_POSITION}")
print(f"N_TEST_PER_CONTEXT: {N_TEST_PER_CONTEXT}")
print(f"N_EXAMPLES_PER_TEST: {N_EXAMPLES_PER_TEST}")
print(f"PSEUDO_COUNTER_2D: {PSEUDO_COUNTER_2D}")
print(f"PSEUDO_COUNTER_3D: {args.pseudo_counter_3d}")


list_nb_example = [N_EXAMPLES_PER_TEST]*N_TEST_PER_CONTEXT

# FOLDER MANAGMENT: global
range_pid_experiment = f"{DATA_RESULT}/Experiment_Brier_PC_3D_SELECTION_{args.pid_inf}_{args.pid_sup}"
list_path_folder = [DATA_RESULT, range_pid_experiment]
for path_folder in list_path_folder:
    if not os.path.exists(path_folder):
        os.makedirs(path_folder,  exist_ok=True)


# MSA PATH
path_folder_seed = f"{DATA_MNHN}/{NAME_FASTA_TEST_FOLDER}"

# PATH 2D_PROBA WITH PSEUDO_COUTER_2D
path_table_2d_proba = f"{DATA_2D_PROBA}/{args.pid_inf}_{args.pid_sup}/proba_{PSEUDO_COUNTER_2D}.npy"
# PATH 3D_PROBA WITH PSEUDO_COUNTER_3D
path_folder_table_3d_proba = f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}/{args.pseudo_counter_3d}"

# PATH PRE-PROECESSING DICO
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

    # dans les expériences réalisées avec ce script,
    #  soit les 4 valeurs de context sont nulles, soit une seule est non nulle
    # et on n'étudiera que vers l'origine ou vers la desination (ne pas mélanger les 2)
    # à contraindre ca? ou je serais la seule utilisatrice et donc je sais quels contextes je peux tester?

    test_number = 0
    total_list_unit_score_brier = []
    for nb_exemple_test in list_nb_example:
        test_number += 1
        print("")
        print("_____________")
        print(f"TEST {test_number}/{N_TEST_PER_CONTEXT}")
        dico_example = selection_example.example_number_per_seed(path_file_dico_seed_normal, nb_exemple_test)
        list_example = selection_example.multi_random_example_selection(path_folder_seed, dico_example,
                                                                        path_folder_dico_seq,
                                                                        context_ol, context_or, context_dl, context_dr,
                                                                        ALPHABET)
    
        list_unit_score_brier, list_unit_score_brier_non_contextuel = brier.brier_score(list_example,
                                                               context_ol, context_or, context_dl, context_dr,
                                                               ALPHABET,
                                                               path_folder_table_3d_proba, path_table_2d_proba,
                                                               args.method)


        # TEST CONVERGENCE
        name_csv = f"convergence_{args.pseudo_counter_3d}"
        path_csv = f"{range_pid_experiment}/{name_csv}.csv"
        file_exists = os.path.isfile(path_csv)

        with open (path_csv, 'a', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            headers = ('context', 'exemple', 'brier_unit', 'brier_unit_no_c')

            if not file_exists:
                writer.writerow(headers)

            for score_brier_unit, score_brier_unit_no_context, ex in zip(list_unit_score_brier, list_unit_score_brier_non_contextuel, list_example):
                data = context, ex,  score_brier_unit, score_brier_unit_no_context
                writer.writerow(data)


        # WORKS ONLY OF ONE QUARTER WINDOW IS CONSIDERED
        if args.reference == "o" and context != (0,0,0,0):
            if context_ol != 0:
                relative_position = -context_ol
            if context_or != 0:
                relative_position = context_or
        if args.reference == "d" and context != (0,0,0,0):
            if context_dl != 0:
                relative_position = -context_dl
            if context_dr != 0:
                relative_position = context_dr
        if context == (0,0,0,0):
            relative_position = 0

        name_csv = f"{N_TEST_PER_CONTEXT}_{nb_exemple_test}_{args.method}_{args.reference}_{args.pseudo_counter_3d}"
        path_selection_pc_3d = f"{range_pid_experiment}/{name_csv}.csv"
        file_exists = os.path.isfile(path_selection_pc_3d)

        # REVIEW REGISTRATION
        with open (path_selection_pc_3d, 'a', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            headers = ('pseudo_counter_3d', 'relative_position', 'score_brier', 'score_brier_ss_context', 'nb_example')

            if not file_exists:
                writer.writerow(headers)

            score_brier = np.mean(list_unit_score_brier)
            score_brier_ss_context = np.mean(list_unit_score_brier_non_contextuel)
            nb_example = len(list_unit_score_brier)

            data = args.pseudo_counter_3d, relative_position, score_brier, score_brier_ss_context, nb_example
            writer.writerow(data)


end = time.time()
print("")
print(f"DONE in {round(end-start, 4)} s")

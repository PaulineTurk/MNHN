"""
Brier Score Experiment:
nohup python3 main_brier.py 30 40 uni origine 0.001> test_brier_$$.txt 2>&1
"""

# IMPORTS
import os
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
parser.add_argument("method", help="'uni' OR 'multi'", type=str)
parser.add_argument("reference", help="'origine' OR 'destination'", type=str)
parser.add_argument("weight_2d", help="weighting 3D with 2D", type=float)
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
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/5_BRIER_SCORE_test_check_method"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
NAME_FASTA_TEST_FOLDER = "Pfam_split/Pfam_test"
MAX_RELATIVE_POSITION = 3
N_TEST_PER_CONTEXT = 1
N_EXAMPLES_PER_TEST = 1_000
PSEUDO_COUNTER_2D = 1
PSEUDO_COUNTER_3D = 0


print("_______________________________________________________________________")
print("                  EVALUATION OF THE PREDICTIVE ABILITY                 ")
print("                     OF A CONTEXTUAL AMINO-ACID                        ")
print("_______________________________________________________________________")

  
print("_________________________________")
print("  Preparation of the experiment  ")
print("_________________________________")

if args.reference == "origine":
    list_ol          = [(i, 0, 0, 0) for i in range(1,MAX_RELATIVE_POSITION + 1)[::-1]]
    list_non_context = [(0, 0, 0, 0)]
    list_or          = [(0, i, 0, 0) for i in range(1,MAX_RELATIVE_POSITION + 1)]

    list_context = list_ol  + list_non_context + list_or


if args.reference == "destination":
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
print(f"PSEUDO_COUNTER_3D: {PSEUDO_COUNTER_3D}")
print(f"WEIGHT_2D: {args.weight_2d}")

list_nb_example = [N_EXAMPLES_PER_TEST]*N_TEST_PER_CONTEXT

# FOLDER MANAGMENT: global
range_pid_experiment = f"{DATA_RESULT}/Experiment_Brier_{args.pid_inf}_{args.pid_sup}"
path_new_file_experiment = f"{range_pid_experiment}/{args.reference}_{args.weight_2d}"
list_path_folder = [DATA_RESULT, range_pid_experiment, path_new_file_experiment]
for path_folder in list_path_folder:
    IS_EXIST = os.path.exists(path_folder )
    if not IS_EXIST:
        os.makedirs(path_folder,  exist_ok=True)


# MSA PATH
path_folder_seed = f"{DATA_MNHN}/{NAME_FASTA_TEST_FOLDER}"

# PATH 2D_PROBA
path_table_2d_proba = f"{DATA_2D_PROBA}/{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER_2D}/proba.npy"
# PATH 3D_PROBA
path_folder_table_3d_proba = f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER_3D}_{args.weight_2d}"

# PATH PRE-PROECESSING DICO
path_folder_dico_seq = f"{DATA_EXEMPLE_TEST}/{args.pid_inf}_{args.pid_sup}/seq"
path_file_dico_seed_normal = f"{DATA_EXEMPLE_TEST}/{args.pid_inf}_{args.pid_sup}/seed/seed_normal.npy"

# INITIALIZATIONS
list_score_brier  = []
list_position     = []

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
        dico_example = selection_example.example_number_per_seed(path_file_dico_seed_normal, N_EXAMPLES_PER_TEST)
        list_example = selection_example.multi_random_example_selection(path_folder_seed, dico_example,
                                                                        path_folder_dico_seq,
                                                                        context_ol, context_or, context_dl, context_dr,
                                                                        ALPHABET)
    
        score_brier, list_unit_score_brier, nb_example, prediction_info = brier.brier_score(list_example,
                                                               context_ol, context_or, context_dl, context_dr,
                                                               ALPHABET,
                                                               path_folder_table_3d_proba, path_table_2d_proba,
                                                               args.method)

        # WORKS ONLY OF ONE QUARTER WINDOW IS CONSIDERE
        if args.reference == "origine" and context != (0,0,0,0):
            if context_ol != 0:
                relative_position = -context_ol
            if context_or != 0:
                relative_position = context_or  
        if args.reference == "destination" and context != (0,0,0,0):
            if context_dl != 0:
                relative_position = -context_dl
            if context_dr != 0:
                relative_position = context_dr
        if context == (0,0,0,0):
            relative_position = 0


        # REVIEW REGISTRATION
        perc_conservation_true, perc_conservation_false, perc_substitution_true, perc_substitution_false = prediction_info
        data = (args.reference, relative_position, test_number, score_brier, nb_example,
                perc_conservation_true, perc_conservation_false,
                perc_substitution_true, perc_substitution_false)

        with open(f"{path_new_file_experiment}/{args.method}.csv", 'a', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            # header: pid_inf, pid_sup, n_valid_seed, n_non_info_seed, n_valid_aa_couple_global,
            #         perc_conservation_true, perc_conservation_false,
            #         perc_substitution_true, perc_substitution_false
            writer.writerow(data)

end = time.time()
diff = round(end-start, 4)
print("")
print(f"DONE in {diff}s")

"""
Brier Score Experiment:
"""

# IMPORTS
import argparse
import numpy as np


import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import brierNeighbour.selection_example as selection_example
import brierNeighbour.



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
MAX_RELATIVE_POSITION = 1
N_TEST_PER_CONTEXT = 1
N_EXAMPLES_PER_TEST = 1
PSEUDO_COUNTER_2D = pow(10, -2)



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

# print(f"PID_INF: {args.pid_inf}")
# print(f"PID_SUP: {args.pid_sup}")
# print(f"METHOD: {args.method}")
# print(f"LIST CONTEXT TESTED: {list_context}")
# print(f"MAX_RELATIVE_POSITION: {MAX_RELATIVE_POSITION}")
# print(f"N_TEST_PER_CONTEXT: {N_TEST_PER_CONTEXT}")
# print(f"N_EXAMPLES_PER_TEST: {N_EXAMPLES_PER_TEST}")
# print(f"PSEUDO_COUNTER_2D: {PSEUDO_COUNTER_2D}")
print(f"PSEUDO_COUNTER_3D: {args.pseudo_counter_3d}")


list_nb_example = [N_EXAMPLES_PER_TEST]*N_TEST_PER_CONTEXT


# MSA PATH
path_folder_seed = f"{DATA_MNHN}/{NAME_FASTA_TEST_FOLDER}"

# PATH 2D_PROBA WITH PSEUDO_COUTER_2D
path_table_2d_proba = f"{DATA_2D_PROBA}/{args.pid_inf}_{args.pid_sup}/proba_{PSEUDO_COUNTER_2D}.npy"
# PATH 3D_PROBA WITH PSEUDO_COUNTER_3D
path_folder_table_3d_proba = f"{DATA_3D_PROBA}/{args.pid_inf}_{args.pid_sup}/{args.pseudo_counter_3d}"

# PATH PRE-PROECESSING DICO
path_folder_dico_seq = f"{DATA_EXEMPLE_TEST}/{args.pid_inf}_{args.pid_sup}/seq"
path_file_dico_seed_normal = f"{DATA_EXEMPLE_TEST}/{args.pid_inf}_{args.pid_sup}/seed/seed_normal.npy"





for context in list_context:
    print(f"\nCONTEXT: {context}")
   
    context_ol, context_or, context_dl, context_dr = context

    test_number = 0
    total_list_unit_score_brier = []
    vector_diff = 0

    for nb_exemple_test in list_nb_example:

        dico_example = selection_example.example_number_per_seed(path_file_dico_seed_normal, N_EXAMPLES_PER_TEST)
        list_example = selection_example.multi_random_example_selection(path_folder_seed, dico_example,
                                                                        path_folder_dico_seq,
                                                                        context_ol, context_or, context_dl, context_dr,
                                                                        ALPHABET)
    
        score_brier, list_unit_score_brier, nb_example, vector_diff, list_cube_quarter_window_ol, list_path_ol= brier.brier_score_check_convergence(list_example,
                                                               context_ol, context_or, context_dl, context_dr,
                                                               ALPHABET,
                                                               path_folder_table_3d_proba, path_table_2d_proba,
                                                               args.method,
                                                               vector_diff)
        




        LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -3),
                          pow(10, -1),
                          pow(10, 0),
                          pow(10, 1),
                          pow(10, 5)]


        table_1d = np.load(f"{DATA_2D_PROBA}/{args.pid_inf}_{args.pid_sup}/freq_1d.npy", allow_pickle='TRUE').item()


        for index, table_3d in enumerate(list_cube_quarter_window_ol):
            print(f"PATH: {list_path_ol[index]}")
            for pseudo_counter_3d in LIST_PSEUDO_COUNTER_3D:
                pseudo_counter_3d = float(pseudo_counter_3d)
                print(f"PSEUDO_COUNTER_3D: {pseudo_counter_3d}")
            
                sum_cases = 0
                for aa_c in ALPHABET:
                    for aa_1 in ALPHABET:
                        for aa_2 in ALPHABET:
                            diff = abs(table_3d[aa_1][aa_2][aa_c] - table_1d[aa_c])
                            #print(f"{aa_1},{aa_2}: {diff}")
                            sum_cases += diff
                print(sum_cases)






    print(f"DIFF: {vector_diff/nb_example}")

# A COMPLETER .............................




# IMPORTS
import os
import argparse
import csv

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import brierNeighbour.selection_example as selection_example
import brierNeighbour.brier as brier
import utils.folder as folder



# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("method", help="'uni' OR 'multi'", type=str)
parser.add_argument("reference", help="'origine' OR 'destination'", type=str)
parser.add_argument("orientation", help="'left' ou 'right'", type=str)
args = parser.parse_args()

# PATH FOR THE MSA ACCESS
DATA_MNHN = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
# PATH FOR THE PRE-PROCESSED DATA TEST EXAMPLES
DATA_EXEMPLE_TEST = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST"
# PATH FOR 2D_PROBA
PATH_FOLDER_2D_PROBA = f"{file.parents[2]}"  # à compléter ................................
# PATH FOR 3D_PROBA
PATH_FOLDER_3D_PROBA = f"{file.parents[2]}"  # à compléter ................................
# PATH FOR THE RESULTS
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/5_BRIER_SCORE"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
NAME_FASTA_TEST_FOLDER = "Pfam_split/Pfam_test"
MAX_RELATIVE_POSITION = 3
N_TEST_EXAMPLES = 10
PSEUDO_COUNTER_2D = 1
WEIGHT_2D = 0.001

print("_______________________________________________________________________")
print("                  EVALUATION OF THE PREDICTIVE ABILITY                 ")
print("                     OF A CONTEXTUAL AMINO-ACID                        ")
print("_______________________________________________________________________")

  
print("_________________________________")
print("  Preparation of the experiment  ")
print("_________________________________")

DICO_VOISINAGE = {("origine", "left"): (MAX_RELATIVE_POSITION, 0, 0, 0),
                  ("origine", "right"): (0, MAX_RELATIVE_POSITION, 0, 0),
                  ("destination", "left"): (0, 0, MAX_RELATIVE_POSITION, 0),
                  ("destination", "right"): (0, 0, 0, MAX_RELATIVE_POSITION)}


context = DICO_VOISINAGE[(args.reference, args.orientation)]
context_ol, context_or, context_dl, context_dr = context
print(f"PID_INF: {args.pid_inf}")
print(f"PID_SUP: {args.pid_sup}")
print(f"CONTEXT: {context}")
print(f"MAX_RELATIVE_POSITION: {MAX_RELATIVE_POSITION}")
print(f"N_TEST_EXAMPLES: {N_TEST_EXAMPLES}")
print(f"PATH_FOLDER_2D_PROBA: {PATH_FOLDER_2D_PROBA}")
print(f"PSEUDO_COUNTER_2D: {PSEUDO_COUNTER_2D}")
print(f"PATH_FOLDER_3D_PROBA: {PATH_FOLDER_3D_PROBA}")
print(f"WEIGHT_2D: {WEIGHT_2D}")


# FOLDER MANAGMENT
new_file_dico = f"Experience_Brier_{args.pid_inf}_{args.pid_sup}"
path_res_folder = f"{DATA_RESULT}/{new_file_dico}"
IS_EXIST = os.path.exists(path_res_folder)
if not IS_EXIST:
    os.makedirs(path_res_folder,  exist_ok=True)


# chemin vers les alignements multiples
path_folder_seed = f"{DATA_MNHN}/{NAME_FASTA_TEST_FOLDER}"


# chemin des dico d'info sur les seed et seq de pré-traitement des exemples
path_folder_dico_seq = f"{DATA_EXEMPLE_TEST}/Exemple_{args.pid_inf}_{args.pid_sup}/seq"   # chemin ou on veut enregistrer les dico_seq
path_file_dico_seed = f"{DATA_EXEMPLE_TEST}/Exemple_{args.pid_inf}_{args.pid_sup}/seed" # chemin ou on veut enregistrer le dico_seed
path_dico_seed_normalised = f"{path_file_dico_seed}/seed_normalised.npy"  # fraction des ex test a prendre par seed








CODE = {"A":0, "R":1, "N":2, "D":3, "C":4,
        "Q":5, "E":6, "G":7, "H":8, "I":9,
        "L":10, "K":11, "M":12, "F":13, "P":14,
        "S":15, "T":16, "W":17, "Y":18, "V":19,
        "B":20,"Z":21, "X":22, "O":23, "U":24,
        "-": 25,
        "*": 26}

def encoding(example, code):
    """Encode a sequence according to the code transformation

    Args:
        example (list): an example of aligned amino-acid with some
                        neighbor amino-acid
        code (dico): dictionary to code each letter in the sequence
                    into a numerical value

    Returns:
        list: a list of floats i.e the numeric equivalent of the initial
              sequence
    """
    list_num = []
    for element in example:
        for char in element:
            list_num.append(int(code[char]))
    return list_num





################################################################################
#                                 Expérience                                   #
################################################################################


dico_fonction_brier = {"uni": brier.brier_score_uni,
                       "multi": brier.brier_score_multi}

print("")
print("##########################")
print("context :", context)
print("##########################")


NAME_CSV = f"data_{args.pid_inf}_{args.pid_sup}_{args.methode}_{args.reference}_{args.sens}"
with open(f"{NAME_CSV}.csv",
    'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)



    dico_example = selection_example.example_number_per_seed(path_dico_seed_normalised, N_TEST_EXAMPLES)
    print("nb ex test demandés:", '{:_}'.format(N_TEST_EXAMPLES))
    list_example = selection_example.multi_random_example_selection(path_folder_seed, dico_example,
                                                                    path_folder_dico_seq,
                                                                    context_kl, context_kr, context_pl, context_pr,
                                                                    ALPHABET)
    for example in list_example:
        print(example)
        data = encoding(example, CODE)
        if data[0] <= 20 and data[1] <= 20:
            print(data)
            writer.writerow(data)

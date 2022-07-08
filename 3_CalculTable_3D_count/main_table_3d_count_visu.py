"""
CONTEXTUAL INFORMATION VIEW
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.autolayout"] = True   # à enlever ?
import argparse
import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0])


# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("ref", help="'o' for origine OR 'd' for destination", type=str)
parser.add_argument("eval",
                    help="'min' for min_count,'not_eval' for not_evaluated OR 'total for total_count",
                    type=str)
args = parser.parse_args()

DATA_RESULT_2D = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
DATA_RESULT_3D = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_COUNT"


ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

# LIST_PID_RANGE = [(0,10), (10,20), (20,30), (30,40), (40,50),
#                   (50,60),(60,70), (70,80), (80,90), (90,100)]
LIST_PID_RANGE = [(0,10)]
INDEX_MAX  = 10

PSEUDO_COUNTER = 1


# FUNCTIONS
def descriptionTable(table_count, alphabet, pseudo_counter, table_type):
    """
    table_count: dico of the table to describe
    table_type: '2D' or '3D'

    return:
    sum over table_count
    number of cells not evaluated in table_count
    median of examples to each cell evaluation
    """
    list_count = []
    total_count = 0
    not_eval_count = 0

    if table_type == "2D":
        for aa_o in alphabet:
            for aa_d in alphabet:
                list_count.append(table_count[aa_o][aa_d])

    if table_type == "3D":
        for aa_o in alphabet:
            for aa_d in alphabet:
                for aa_c in alphabet:
                    list_count.append(table_count[aa_o][aa_d][aa_c])

    for count in list_count:
        total_count += count
        if count == pseudo_counter:
            not_eval_count += 1

    not_eval_percent = 100*not_eval_count/(len(alphabet)**3)
    min_count = min(list_count)

    return total_count, not_eval_percent, min_count


DICO_TITLE_X = {"o": "Position de l'aa contextuel par rapport à l'aa d'origine",
                "d": "Position de l'aa contextuel par rapport à l'aa de destination"}

DICO_TITLE_Y = {"not_eval": "Pourcentage de paramètres non-estimés",
                "total": "Nombre d'exemples d'apprentissage",
                "min":"Nombre min d'exemples d'apprentissage par paramètre"}
                
DICO_TITLE_GRAPH = {"not_eval":"Pourcentage de paramètres non estimés \ndans les tables de comptage",
                    "total": "Nombre des exemples d'apprentissage \npour le calcul des tables de comptage",
                    "min": "Nombre minimal d'exemples d'apprentissage par paramètre\n pour le calcul des tables de comptage"}

plt.figure(figsize=[10,9])

list_position = [i for i in range(-INDEX_MAX,INDEX_MAX+1)]
for pid_range in LIST_PID_RANGE:
    pid_inf, pid_sup = pid_range
    total_count_global, not_eval_percent_global, min_count_global = [], [], []

    DICO_ESTIMATION = {"not_eval": not_eval_percent_global,
                        "total": total_count_global,
                        "min": min_count_global}

    # chemin des tables 3d count par tranche de pid
    path_table_3d = f"{DATA_RESULT_3D}/table_3d_count_{pid_inf}_{pid_sup}_{PSEUDO_COUNTER}"

    
    # importation des tables de comptage 3D de la tranche en position négative
    list_neg_index = [i for i in range(1,INDEX_MAX+1)][::-1]
    for index in list_neg_index:
        table3d_comptage = np.load(f"{path_table_3d}/table_3d_count_({-index},{args.ref}).npy", allow_pickle='TRUE').item()
        total_count, not_eval_percent, min_count = descriptionTable(table3d_comptage, ALPHABET,
                                                                    PSEUDO_COUNTER, "3D")
        total_count_global.append(total_count)
        not_eval_percent_global.append(not_eval_percent)
        min_count_global.append(min_count)

    # importation de la table de comptage 2D de la tranche
    path_table_2d    = f"{DATA_RESULT_2D}/table_2d_{pid_inf}_{pid_sup}_{PSEUDO_COUNTER}"
    table2d_comptage = np.load(f"{path_table_2d}/table_2d_count.npy", allow_pickle='TRUE').item()
    total_count, not_eval_percent, min_count = descriptionTable(table2d_comptage, ALPHABET,
                                                                PSEUDO_COUNTER, "2D")
    total_count_global.append(total_count)
    not_eval_percent_global.append(not_eval_percent)
    min_count_global.append(min_count)

    # importation des tables de comptage 3D de la tranche en position positive
    for index in range(1,INDEX_MAX+1):
        table3d_comptage = np.load(f"{path_table_3d}/table_3d_count_({index},{args.ref}).npy", allow_pickle='TRUE').item()
        total_count, not_eval_percent, min_count = descriptionTable(table3d_comptage, ALPHABET,
                                                                    PSEUDO_COUNTER, "3D")
        total_count_global.append(total_count)
        not_eval_percent_global.append(not_eval_percent)
        min_count_global.append(min_count)



    plt.plot(list_position, DICO_ESTIMATION[args.eval], label = f"{pid_inf}, {pid_sup}")
    plt.xticks(range(-INDEX_MAX,INDEX_MAX+1))


plt.xlabel(DICO_TITLE_X[args.ref], fontsize=13)
plt.ylabel(DICO_TITLE_Y[args.eval], fontsize=13)
if args.eval == "not_eval":
    plt.yticks(np.arange(0,60,5))
else:
    plt.yscale('log')
plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)                                        
plt.title(DICO_TITLE_GRAPH[args.eval], loc='center',  fontsize=20)
plt.legend()
plt.savefig(f"{args.eval}_{args.ref}.png")

print("DONE")

"""
VIEW: BRIER TEST - MULTI
"""

# IMPORTS

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


file = Path(__file__).resolve()
sys.path.append(file.parents[0])


# PARAMETERS
# EXPERIMENT_NAME = "EXP_1M_UNI"
EXPERIMENT_NAME = "EXP_1M_MULTI"
DATA = f"{file.parents[2]}/MNHN_RESULT/6_TEST_BRIER_NAIVE_BAYES/{EXPERIMENT_NAME}"

L = 6
origin_destination = "origine"
# origin_destination = "destination"
# uni_multi = "uni"
uni_multi = "multi"

LIST_PSEUDO_COUNTER_3D = [pow(10, -2)]

# PROGRAM

# LOAD RESULT
dict_dict_3D = {}
for direction in ["ol", "or", "dl", "dr"]:
    dict_dict_3D[direction] = np.load(f"{DATA}/SCORE_{direction}_{uni_multi}.npy", allow_pickle='True').item()


# GET X, Y

plt.figure(figsize=(8, 6), dpi=80)

list_position_origin = []
for pseudo_counter_3D in LIST_PSEUDO_COUNTER_3D:
    list_position_origin = []
    list_score_origin = []
    if origin_destination == "origine":
        for direction in ["ol", "or"]:
            if direction == "ol":
                for position in range(L, -1, -1):
                    list_position_origin.append(-position)
                    list_score_origin.append(np.mean(dict_dict_3D[direction][position][pseudo_counter_3D]))
            if direction == "or":
                for position in range(0, L+1):
                    list_position_origin.append(position)
                    list_score_origin.append(np.mean(dict_dict_3D[direction][position][pseudo_counter_3D]))
    if origin_destination == "destination":
        for direction in ["dl", "dr"]:
            if direction == "dl":
                for position in range(L, -1, -1):
                    list_position_origin.append(-position)
                    list_score_origin.append(np.mean(dict_dict_3D[direction][position][pseudo_counter_3D]))
            if direction == "dr":
                for position in range(0, L+1):
                    list_position_origin.append(position)
                    list_score_origin.append(np.mean(dict_dict_3D[direction][position][pseudo_counter_3D]))
            
    plt.plot(list_position_origin, list_score_origin, label=pseudo_counter_3D)

plt.legend()
plt.xlabel("RELATIVE POSITION", fontsize=13)
plt.xticks(range(- L, L + 1))
plt.ylabel('BRIER SCORE', fontsize=13)
plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
title_fig = f"Calcul du score de Brier avec bayésien naif sur Pfam_TEST \nréférence: {origin_destination}, méthode: {uni_multi}"
plt.title(title_fig , loc='center', fontsize=20)
title_fig = f"PC_3D_{origin_destination}_{uni_multi}"
plt.savefig(f"{DATA}/{title_fig}.png")

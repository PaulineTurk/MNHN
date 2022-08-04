"""
VIEW: SELECTION PSEUDO_COUNTER 3D ORIGIN & DESTINATION TEST
"""

# IMPORTS

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import itertools


file = Path(__file__).resolve()
sys.path.append(file.parents[0])


# PARAMETERS
# EXPERIMENT_NAME = "EXP_1M_UNI"
EXPERIMENT_NAME = "EXP_1M_MULTI"
DATA = f"{file.parents[2]}/MNHN_RESULT/5_PC_3D_SELECTION_TEST/{EXPERIMENT_NAME}"

L = 6


# uni_multi = "uni"
uni_multi = "multi"


LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10,  0),
                          pow(10,  1),
                          pow(10,  2),
                          pow(10,  3)]

# PROGRAM
import matplotlib as mpl
def set_color_cycle(self, clist=None):
    if clist is None:
        clist = mpl.rcParams['axes.color_cycle']
    self.color_cycle = itertools.cycle(clist)

# LOAD RESULT
dict_dict_3D = {}
for direction in ["ol", "or", "dl", "dr"]:
    dict_dict_3D[direction] = np.load(f"{DATA}/SCORE_{direction}_{uni_multi}.npy", allow_pickle='True').item()

plt.figure(figsize=(10, 8), dpi=100)

for origin_destination in ["origine", "destination"]:
    # GET X, Y
    plt.gca().set_prop_cycle(None)



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
        if origin_destination == "origine":
            # line1, = plt.plot(list_position_origin, list_score_origin, label=pseudo_counter_3D)
            plt.plot(list_position_origin, list_score_origin, label=pseudo_counter_3D)
        if origin_destination == "destination":
            # line2, = plt.plot(list_position_origin, list_score_origin, label=pseudo_counter_3D, linestyle='dashed')
            plt.plot(list_position_origin, list_score_origin, linestyle='dashed')


plt.legend(loc='lower right')
plt.xlabel("RELATIVE POSITION", fontsize=13)
plt.xticks(range(- L, L + 1))
plt.ylabel('BRIER SCORE', fontsize=13)
plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
if uni_multi == "uni":
    title_fig = f"Calcul du score de Brier bayésien naïf avec pseudo-compte 3D \n1M d'exemples TEST, méthode: UNI"
if uni_multi == "multi":
    title_fig = f"Calcul du score de Brier bayésien naïf avec pseudo-compte 3D \n1M d'exemples TEST, méthode: MULTI"
plt.title(title_fig , loc='center', fontsize=20)
title_fig = f"PC_3D_{uni_multi}"
plt.savefig(f"{DATA}/{title_fig}.png")



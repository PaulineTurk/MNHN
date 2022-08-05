"""
VIEW: BRIER TEST - FULL
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
EXPERIMENT_NAME = "EXP_10M_FULL"
DATA = f"{file.parents[2]}/MNHN_RESULT/6_TEST_BRIER_NAIVE_BAYES/{EXPERIMENT_NAME}"

L = 6

LIST_PSEUDO_COUNTER_3D = [pow(10, 0)]

# PROGRAM

# LOAD RESULT
dict_3D = np.load(f"{DATA}/SCORE_FULL.npy", allow_pickle='True').item()


# GET X, Y

plt.figure(figsize=(9, 8), dpi=400)

list_position_origin = []
for pseudo_counter_3D in LIST_PSEUDO_COUNTER_3D:
    list_position = []
    list_score = []

    for position in range(0, L+1):
        list_position.append(position)
        list_score.append(np.mean(dict_3D[position][pseudo_counter_3D]))
          
            
    plt.plot(list_position, list_score, label=pseudo_counter_3D)

# plt.legend()
plt.xlabel("POSITION RELATIVE", fontsize=13)
plt.xticks(range(0, L + 1))
plt.ylabel('SCORE DE BRIER', fontsize=13)
plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
title_fig = f"Calcul du score de Brier bayésien naïf sur Pfam_TEST \nméthode FULL"
plt.title(title_fig , loc='center', fontsize=20)
title_fig = f"PC_3D_FULL_10M"
plt.savefig(f"{DATA}/{title_fig}.png")

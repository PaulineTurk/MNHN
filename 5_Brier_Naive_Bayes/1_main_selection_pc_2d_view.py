"""
VIEW: NON-CONTEXTUAL BRIER SCORE EXPERIMENT
"""

# IMPORTS

import os
import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


file = Path(__file__).resolve()
sys.path.append(file.parents[0])



# PARAMETERS
EXPERIMENT_NAME = "EXP_1M"
DATA = f"{file.parents[2]}/MNHN_RESULT/5_PC_2D_SELECTION/{EXPERIMENT_NAME}"


# PROGRAM

# LOAD RESULT
dict_score = np.load(f"{DATA}/SCORE.npy", allow_pickle='True').item()

# GET THE PSEUDO-COUNTER 2D VALUES
list_pseudo_counter_2D = dict_score.keys()

# GET THE BRIER SCORES
list_mean_score_Brier = []
for pseudo_counter_2D in list_pseudo_counter_2D:
    list_mean_score_Brier.append(np.mean(dict_score[pseudo_counter_2D]))


# VIEW
plt.plot(list_pseudo_counter_2D, list_mean_score_Brier)
plt.xlabel("PSEUDO-COUNTER 2D", fontsize=13)
plt.xscale('log')
plt.ylabel('BRIER SCORE', fontsize=13)
plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
title_fig = "Selection du pseudo-compte des tables_2d"
plt.title(title_fig , loc='center', fontsize=20)
title_fig = "PC_2D"
plt.savefig(f"{DATA}/{title_fig}.png")

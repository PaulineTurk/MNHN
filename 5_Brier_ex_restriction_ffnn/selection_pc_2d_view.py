"""
Visualisation of the Brier Score tests:
/home/pauline/Bureau/MNHN_RESULT/5_BRIER_SCORE/Experiment_Brier_60_70/origine_0.1
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
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
args = parser.parse_args()



DATA = f"{file.parents[2]}/MNHN_RESULT/5_PC_2D_SELECTION/Experiment_Brier_PC_2D_SELECTION"



print("_______________________________________________________________________")
print("                      VIEW PSEUDO_COUNTER_2D                           ")
print("_______________________________________________________________________")

# FOLDER REGISTRATION
path_image = f"{DATA}/REVIEW"
if not os.path.exists(path_image):
    os.makedirs(path_image)

# CSV READING IN DATAFRAME
# ADD FOR LOOP ON PID RANGE
path_csv_file = f"{DATA}/{args.pid_inf}_{args.pid_sup}.csv"
# conversion to pandas dataframe + headers addition
df = pd.read_csv(path_csv_file,
                 names=['pseudo_counter_2d', 'score_brier', 'nb_example'])

df_mean = df.groupby('pseudo_counter_2d')['score_brier'].mean()
print(df_mean)
list_pseudo_counter_2d = np.array([float(elem) for elem in df_mean.keys()])
print(list_pseudo_counter_2d)
list_mean = np.array([elem for elem in df_mean])
df_std = df.groupby('pseudo_counter_2d')['score_brier'].std()
list_std = np.array([elem for elem in df_std])


plt.plot(list_pseudo_counter_2d, list_mean)
plt.fill_between(list_pseudo_counter_2d, list_mean - 2*list_std,
                     list_mean + 2*list_std,
                     alpha=0.5, label = "pseudo_counter_2d")

plt.xlabel("PSEUDO-COUNTER 2D", fontsize=13)
plt.xscale('log')
plt.ylabel('BRIER SCORE', fontsize=13)
plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
title = "Selection du pseudo-compte des tables_2d"
plt.title(title, loc='center', fontsize=20)
plt.legend(f"{args.pid_inf}_{args.pid_sup}")
title = "PC_2D"
plt.savefig(f"{path_image}/{title}.png")

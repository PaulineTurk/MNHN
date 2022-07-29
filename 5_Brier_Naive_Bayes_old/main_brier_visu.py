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
from matplotlib import rc

# latex activation
rc('text', usetex=True)
# sudo apt-get install dvipng texlive-latex-extra texlive-fonts-recommended  cm-super

file = Path(__file__).resolve()
sys.path.append(file.parents[0])




# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("method", help="'uni' OR 'multi'", type=str)
parser.add_argument("reference", help="'origine' OR 'destination'", type=str)
args = parser.parse_args()



DATA = f"{file.parents[2]}/MNHN_RESULT/5_BRIER_SCORE/Experiment_Brier_{args.pid_inf}_{args.pid_sup}"
LIST_WEIGHT = [0.001, 0.01, 0.1, 1.0, 10.0]

print("_______________________________________________________________________")
print("                           WEIGHT 2D EFFECT                            ")
print("_______________________________________________________________________")

# FOLDER REGISTRATION
path_image = f"{DATA}/REVIEW"
IS_EXIST = os.path.exists(path_image)
if not IS_EXIST:
    os.makedirs(path_image)

# CSV READING IN DATAFRAME
# for weight in LIST_WEIGHT:
for weight in LIST_WEIGHT:
    relative_path_csv_file = f"{args.reference}_{weight}/{args.method}.csv"
    # conversion to pandas dataframe + headers addition
    df = pd.read_csv(f"{DATA}/{relative_path_csv_file}",
                 names=['o_d', 'relative_position', 'test_num', 'Brier_Score', 'n_examples',
                         'perc_conservation_true', 'perc_conservation_false',
                         'perc_substitution_true', 'perc_substitution_false'])

    df_mean = df.groupby('relative_position')['Brier_Score'].mean()
    list_relative_postition = np.array([elem for elem in df_mean.keys()])
    list_mean = np.array([elem for elem in df_mean])
    df_std = df.groupby('relative_position')['Brier_Score'].std()
    list_std = np.array([elem for elem in df_std])


    plt.plot(list_relative_postition, list_mean)
    plt.fill_between(list_relative_postition, list_mean - 2*list_std,
                     list_mean + 2*list_std,
                     alpha=0.5, label = f"{weight}")
    max_position = max(list_relative_postition)
    plt.xticks(range(- max_position, max_position + 1))
    #plt.yticks(np.arange(0,1.1,0.1))

    plt.xlabel("Position contextuelle relative", fontsize=13)
    plt.ylabel('Score de Brier', fontsize=13)
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
    plt.title(f"{args.pid_inf}_{args.pid_sup}_{args.method}_{args.reference}",
                loc='center', fontsize=20)
    plt.legend(title=r"\textbf{POIDS 2D}")
    plt.savefig(f"{path_image}/{args.pid_inf}_{args.pid_sup}_{args.method}_{args.reference}.png")

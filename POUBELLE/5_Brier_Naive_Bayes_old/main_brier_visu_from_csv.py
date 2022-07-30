"""
Visualisation of the Brier Score tests:
python3 main_brier_visu_from_csv.py 60 70 uni o
"""

# IMPORTS



import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


file = Path(__file__).resolve()
sys.path.append(file.parents[0])


# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
parser.add_argument("method", help="'uni' OR 'multi'", type=str)
parser.add_argument("reference", help="'origine' OR 'destination'", type=str)
args = parser.parse_args()

EXP_NAME = "5_PC_3D_SELECTION_CORRECTION_BIAIS_test_3"
MAIN_DATA = f"{file.parents[2]}/MNHN_RESULT/{EXP_NAME}/Experiment_Brier_PC_3D_SELECTION_{args.pid_inf}_{args.pid_sup}"

N_TEST_PER_CONTEXT = 1
N_EX_PER_TEST = 1_000_000


LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -5),
                          pow(10, -4),
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10, 0),
                          pow(10, 1),
                          pow(10, 2),
                          pow(10, 3),
                          pow(10, 4),
                          pow(10, 5),
                          pow(10, 6),
                          pow(10, 7),
                          pow(10, 8),
                          pow(10, 9),
                          pow(10, 10)]



# FOLDER REGISTRATION
path_image = f"{MAIN_DATA}/REVIEW_ON_FULL_SCORE"
if not os.path.exists(path_image):
    os.makedirs(path_image)

f,a = plt.subplots(3,6)
a = a.ravel()
index = 0
for ax in a:
    if index < len(LIST_PSEUDO_COUNTER_3D):
        PC = float(LIST_PSEUDO_COUNTER_3D[index])

        DATA = f"{MAIN_DATA}/{N_TEST_PER_CONTEXT}_{N_EX_PER_TEST}_{args.method}_{args.reference}_{PC}.csv"

        df = pd.read_csv(DATA,
        names=['pseudo_counter_3d', 'relative_position', 'score_brier', 'score_brier_ss_context', 'nb_example'])
        
        x = [int(elem) for elem in df['relative_position'][1:].tolist()]
        y_c = [float(elem) for elem in df['score_brier'][1:].tolist()]
        y_no_c = [float(elem) for elem in df['score_brier_ss_context'][1:].tolist()]

        bins = np.linspace(0, 2, 100)


        ax.plot(x, y_c, alpha=0.5, label='C')
        ax.plot(x, y_no_c, alpha=0.5, label='/C')
        ax.legend(loc='upper right')
        ax.set_title(f"{PC}")
        # ax.set_title(titles[idx])
        # ax.set_xlabel(xaxes[idx])
        # ax.set_ylabel(yaxes[idx])
        index += 1
# plt.tight_layout()

# deleate the last empty graph
plt.xlabel("Position")
plt.ylabel("Brier Score")
f.delaxes(ax)

plt.show()


# plt.savefig(f"{path_image}/{title_fig}.png")

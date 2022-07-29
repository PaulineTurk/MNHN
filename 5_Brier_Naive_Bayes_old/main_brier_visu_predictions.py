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
from plotly import graph_objects as go

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



DATA = f"{file.parents[2]}/MNHN_RESULT/5_BRIER_SCORE_test_check_method/Experiment_Brier_{args.pid_inf}_{args.pid_sup}"
LIST_WEIGHT = [0.001]

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
    print(df)

    df_mean_n_examples = df.groupby('relative_position')['n_examples'].mean()
    list_relative_postition = np.array([int(elem) for elem in df_mean_n_examples.keys()])
    list_mean_n_examples = [int(elem) for elem in df_mean_n_examples]
    names_x = []
    for position, n_ex in zip(list_relative_postition, list_mean_n_examples):
        names_x.append(f"{position};{n_ex}")


    df_mean_conservation_true = df.groupby('relative_position')['perc_conservation_true'].mean()
    list_mean_conservation_true = np.array([elem for elem in df_mean_conservation_true])

    df_mean_conservation_false = df.groupby('relative_position')['perc_conservation_false'].mean()
    list_mean_conservation_false = np.array([elem for elem in df_mean_conservation_false])

    df_mean_substitution_true = df.groupby('relative_position')['perc_substitution_true'].mean()
    list_mean_substitution_true = np.array([elem for elem in df_mean_substitution_true])

    df_mean_substitution_false = df.groupby('relative_position')['perc_substitution_false'].mean()
    list_mean_substitution_false = np.array([elem for elem in df_mean_substitution_false])





fig = go.Figure(
    data=[
        go.Bar(
            name="CONSERVATION_VRAI",
            x=names_x,
            y=list_mean_conservation_true,
            offsetgroup=0,
        ),
        go.Bar(
            name="CONSERVATION_FAUX",
            x=names_x,
            y=list_mean_conservation_false,
            offsetgroup=0,
            base=list_mean_conservation_true,
        ),
        go.Bar(
            name="SUBSTITUTION_VRAI",
            x=names_x,
            y=list_mean_substitution_true,
            offsetgroup=2,
        ),
        go.Bar(
            name="SUBSTITUTION_FAUX",
            x=names_x,
            y=list_mean_substitution_false,
            offsetgroup=2,
            base=list_mean_substitution_true,
        )
    ],
    layout=go.Layout(
        title=f"Pourcentage prédiction {args.pid_inf}_{args.pid_sup}",
        yaxis_title="Pourcentage prédiction"
    )
)

fig.show()

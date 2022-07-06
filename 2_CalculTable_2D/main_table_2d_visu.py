"""
NON-CONTEXTUAL VIEW
"""

# IMPORTS
import sys
import os
import argparse
from pathlib import Path
import numpy as np

file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table2d.table2dfonction as table2dfonction
import utils.folder as folder


# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
args = parser.parse_args()

DATA = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
NAME_FASTA_TRAIN_FOLDER = "Pfam_split/Pfam_train"
NAME_PID_FOLDER         = "PID"
PSEUDO_COUNTER = 1




print("_______________________________________________________________________")
print("                           VIEW 2D_TABLES                              ")
print("                          BY RANGE OF PID                              ")
print("_______________________________________________________________________")

path_folder_fasta      = f"{DATA}/{NAME_FASTA_TRAIN_FOLDER}"
path_folder_pid        = f"{DATA}/{NAME_PID_FOLDER}"
path_res = f"{DATA_RESULT}/table_2d_{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER}"
path_res_graph = folder.creat_folder(f"{path_res}/graphe")

# table_2d score (BLOSUM formula)
print(f"\ntable_2d_score ({args.pid_inf},{args.pid_sup}) :\n")
table_2d_score =  np.load(f"{path_res}/table_2d_score.npy", allow_pickle='TRUE').item()
table2dfonction.table_2d_visualisation_transposition(table_2d_score)
name_folder_fasta =  os.path.basename(path_folder_fasta)
title_heatmap = f"Heatmap de la table_2d de scores [{args.pid_inf},{args.pid_sup}] calculée sur {name_folder_fasta}"
table2dfonction.table_2d_heatmap(table_2d_score, path_res_graph, title_heatmap, size_annot = 5)


# table_2d score comparison with BLOSUM_62
PID_INF_REF = 62
matrix_diff, PID_INF_REF, average_diff = table2dfonction.table_2d_difference(table_2d_score, PID_INF_REF)
title_heatmap = f"Heatmap des différences entre la table_2d de scores [{args.pid_inf},{args.pid_sup}] et la Blosum_{PID_INF_REF} de référence\nLa différence moyenne de score est de : {average_diff}"
table2dfonction.table_2d_heatmap(matrix_diff, path_res_graph, title_heatmap, size_annot = 5)


# table_2d probability
print(f"\ntable_2d_proba ({args.pid_inf},{args.pid_sup}) :\n")
table_2d_proba =  np.load(f"{path_res}/table_2d_proba.npy", allow_pickle='TRUE').item()

table2dfonction.table_2d_visualisation_transposition(table_2d_proba)
table2dfonction.sum_line_transposition(table_2d_proba)
title_heatmap = f"Heatmap de la table_2d de probabilités conditionnelles [{args.pid_inf},{args.pid_sup}] calculée sur {name_folder_fasta}"
table2dfonction.table_2d_heatmap(table_2d_proba, path_res_graph, title_heatmap, size_annot = 3)

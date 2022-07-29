"""
NON-CONTEXTUAL VIEW:
python3 main_table_2d_view.py 90 100 > data_table_2d_view_90_100_$$.txt 2>&1
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


# PARAMETERS
L = 6
PID_INF = 40
PID_INF = 50

DATA = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D"
NAME_FASTA_TRAIN_FOLDER = "Pfam_split/Pfam_train"
NAME_PID_FOLDER = "PID"
ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
LIST_PSEUDO_COUNTER_2D = [0,
                          pow(10, -5),
                          pow(10, -4),
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          1,
                          10]




# PROGRAM

path_folder_fasta      = f"{DATA}/{NAME_FASTA_TRAIN_FOLDER}"
path_folder_pid        = f"{DATA}/{NAME_PID_FOLDER}"
path_res = f"{DATA_RESULT}/{args.pid_inf}_{args.pid_sup}"
path_res_graph = f"{path_res}/graph"
if not os.path.exists(path_res_graph):
    os.makedirs(path_res_graph)

# table_2d score (BLOSUM formula)
print(f"\ntable_2d_score ({args.pid_inf},{args.pid_sup}) :\n")
table_2d_score = np.load(f"{path_res}/score.npy", allow_pickle='TRUE').item()
table2dfonction.table_2d_visualisation(table_2d_score)
name_folder_fasta =  os.path.basename(path_folder_fasta)
title_heatmap = f"Heatmap de la table_2d de scores [{args.pid_inf},{args.pid_sup}] calculée sur {name_folder_fasta}"
table2dfonction.table_2d_heatmap(table_2d_score, path_res_graph, title_heatmap, size_annot = 5)


# table_2d score comparison with BLOSUM_62
PID_INF_REF = 62
matrix_diff, PID_INF_REF, average_diff = table2dfonction.table_2d_difference(
                                                            table_2d_score,
                                                            ALPHABET,
                                                            PID_INF_REF)
title_heatmap = f"Heatmap des différences entre la table_2d de scores [{args.pid_inf},{args.pid_sup}] et la Blosum_{PID_INF_REF} de référence\nLa différence moyenne de score est de : {average_diff}"
table2dfonction.table_2d_heatmap(matrix_diff, path_res_graph, title_heatmap, size_annot = 5)


# table_2d probability
for pseudo_counter_2d in LIST_PSEUDO_COUNTER_2D:
    table_2d_proba =  np.load(f"{path_res}/proba_{pseudo_counter_2d}.npy", allow_pickle='TRUE').item()
    table2dfonction.table_2d_visualisation(table_2d_proba)
    table2dfonction.sum_line(table_2d_proba)
    title_heatmap = f"Heatmap de la table_2d de probabilités conditionnelles [{args.pid_inf},{args.pid_sup}] calculée sur {name_folder_fasta}\n pseudo_counter_2d: {pseudo_counter_2d}"
    table2dfonction.table_2d_heatmap(table_2d_proba, path_res_graph, title_heatmap, size_annot = 3)

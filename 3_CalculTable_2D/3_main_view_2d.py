"""
2D_VIEW
nohup python3 3_main_view_2d.py > 3.out 2>&1 &
"""

# IMPORTS
import sys
import os
from pathlib import Path
import numpy as np
import time
from datetime import datetime

file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table2d.function_table2d as function_table2d


# PARAMETERS
PID_INF = 40
PID_SUP = 50
L = 6

DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

LIST_PSEUDO_COUNTER_2D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10,  0),
                          pow(10,  1),
                          pow(10,  2),
                          pow(10,  3)]




# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

start = time.time()

path_res_graph = f"{DATA_RESULT}/GRAPH"
os.makedirs(path_res_graph, exist_ok=True)

# table_2d score (BLOSUM formula)
print(f"\ntable_2d_score ({PID_INF},{PID_SUP}) :\n")
table_2d_score = np.load(f"{DATA_RESULT}/score.npy", allow_pickle='TRUE').item()
function_table2d.table_2d_visualisation(table_2d_score)
title_heatmap = f"Heatmap de la table_2d de scores [{PID_INF},{PID_SUP}] calculée sur Pfam TRAIN"
function_table2d.table_2d_heatmap(table_2d_score, path_res_graph, title_heatmap, size_annot = 5)


# table_2d score comparison with BLOSUM_62
PID_INF_REF = 62
matrix_diff, PID_INF_REF, average_diff = function_table2d.table_2d_difference(
                                                            table_2d_score,
                                                            ALPHABET,
                                                            PID_INF_REF)
title_heatmap = f"Heatmap des différences entre la table_2d de scores [{PID_INF},{PID_SUP}] et la Blosum_{PID_INF_REF} de référence\nLa différence moyenne de score est de : {average_diff}"
function_table2d.table_2d_heatmap(matrix_diff, path_res_graph, title_heatmap, size_annot = 5)


# table_2d probability
for pseudo_counter_2d in LIST_PSEUDO_COUNTER_2D:
    print("")
    print(f"\ntable_2d_proba [{PID_INF},{PID_SUP}]")
    print(f"PSEUDO_COUNTER_2D: {pseudo_counter_2d}")
    print("")
    pseudo_counter_2d = float(pseudo_counter_2d)
    table_2d_proba =  np.load(f"{DATA_RESULT}/proba_{pseudo_counter_2d}.npy", allow_pickle='TRUE').item()
    function_table2d.table_2d_visualisation(table_2d_proba)
    function_table2d.sum_line(table_2d_proba)
    title_heatmap = f"Heatmap de la table_2d de probabilités conditionnelles [{PID_INF},{PID_SUP}] calculée sur Pfam TRAIN\n pseudo_counter_2d: {pseudo_counter_2d}"
    function_table2d.table_2d_heatmap(table_2d_proba, path_res_graph, title_heatmap, size_annot = 3)

end = time.time()
print(f"\nDONE IN: {round(end - start, 2)} s")

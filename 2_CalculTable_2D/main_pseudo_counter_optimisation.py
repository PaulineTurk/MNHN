"""
PSEUDO_COUNTER OPTIMISATION:
python3 main_pseudo_counter_optimisation.py 90 100 > data_pc_2d_optimisation_$$.txt 2>&1
"""


# IMPORTS
import argparse
import csv
import numpy as np

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])
import table2d.table2dfonction as table2dfonction


# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur", type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur", type=int)
args = parser.parse_args()

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

DATA = f"{file.parents[2]}/MNHN_RESULT_DRAFT/2_TABLE_2D"
PSEUDO_COUNTER_INITIAL = 1
LIST_PSEUDO_COUNTER_2D = [0, 1, 10, 100, 1_000, 10_000, 100_000, 1_000_000]


print("_______________________________________________________________________")
print("              SELECTION OF PSEUDO_COUNTER_2D BASED ON                  ")
print("                      THE MIN OF PROBA ESTIMATED                       ")
print("_______________________________________________________________________")


with open(f"{DATA}/PC_optimisation.csv", 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header = ("pseudo_counter_2d", "min_proba_estimated")
    writer.writerow(header)
    for pseudo_counter_2d in LIST_PSEUDO_COUNTER_2D:
        path_table_2d = f"{DATA}/{args.pid_inf}_{args.pid_sup}_{PSEUDO_COUNTER_INITIAL}"
        table_2d_proba =  np.load(f"{path_table_2d}/proba_{pseudo_counter_2d}.npy", allow_pickle='TRUE').item()
        proba_min = table2dfonction.min_2D(table_2d_proba, ALPHABET)
        data = (pseudo_counter_2d, proba_min)
        writer.writerow(data)



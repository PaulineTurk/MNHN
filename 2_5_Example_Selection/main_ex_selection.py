"""
Preprocessing of data selection:
nohup python3 main_ex_selection.py > main_ex_selection.out 2>&1 &
"""

# IMPORTS
import time
import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import data_exemple_test.function_ex_selection as ex_selection


DIRECTION = "dr"
NAME_ALL_EXAMPLES_CSV = "test_concat"
PID_INF, PID_SUP = 40, 50
L = 6


DATA = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES"
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/EXAMPLES_{L}_{PID_INF}_{PID_SUP}_{DIRECTION}"
path_original = f"{DATA}/{NAME_ALL_EXAMPLES_CSV}.csv"      # all_examples for max position context fixed & on pid_inf, pid_sup fixed
path_target = f"{DATA}/{NAME_ALL_EXAMPLES_CSV}_{DIRECTION}.csv"
os.makedirs(path_target, exist_ok=True)

start = time.time()

if DIRECTION == "ol":
    ex_selection.ex_selection_ol(path_original, path_target, L)
if DIRECTION == "or":
    ex_selection.ex_selection_or(path_original, path_target, L)
if DIRECTION == "dl":
    ex_selection.ex_selection_dl(path_original, path_target, L)
if DIRECTION == "dr":
    ex_selection.ex_selection_dr(path_original, path_target, L)

end = time.time()

print(f"DONE IN: {round(end - start, 2)} s")
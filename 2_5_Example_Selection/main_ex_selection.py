"""
Preprocessing of data selection:
bash main_ex_selection.sh > main_ex_selection.out 2>&1 &
"""

# IMPORTS

import os
import os.path
import argparse

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import data_exemple_test.function_ex_selection as ex_selection
import utils.folder as folder


# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("path_example_csv",
                     help="path of file to describe the seq info from", type=str)
args = parser.parse_args()


DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT_MINI/2_5_EXAMPLES"
original = args.path_example_csv
path_target = f"{DATA_RESULT}/EXAMPLES_6_restricted_ol"
accession_num = folder.get_accession_number(original)
os.makedirs(path_target, exist_ok=True)

target = f'{path_target}/{accession_num}.csv'

pid_inf, pid_sup = 30, 40
L = 6



ex_selection.ex_selection_ol(original, pid_inf, pid_sup, target, L)
# ex_selection.ex_selection_or(original, pid_inf, pid_sup, target, L)
# ex_selection.ex_selection_dl(original, pid_inf, pid_sup, target, L)
# ex_selection.ex_selection_dr(original, pid_inf, pid_sup, target, L)
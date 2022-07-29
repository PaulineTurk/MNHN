"""
Preprocessing of data selection:
bash 2_main_ex_save.sh > 2_main_ex_save.out 2>&1 &
"""

# IMPORTS

import os
import os.path
import argparse
import numpy as np

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import data_exemple_test.function_ex_save as function_ex_save
import utils.folder as folder
import utils.fastaReader as fastaReader

# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("path_fasta_file",
                     help="path of file to describe the seq info from", type=str)
args = parser.parse_args()


DATA = args.path_fasta_file
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/2_EXAMPLES"
DATA_SEQ_INFO = f"{DATA_RESULT}/SEQ_INFO"
DATA_PID = f"{file.parents[2]}/MNHN_RESULT/1_DATA/PID"

L = 6
PID_INF = 40
PID_SUP = 50

# FOLDER MANAGEMENT
new_folder_example = f"{DATA_RESULT}/EXAMPLES_{L}_{PID_INF}_{PID_SUP}"

list_folder = [new_folder_example]
for path_folder in list_folder:
    os.makedirs(path_folder, exist_ok=True)


# SELECTION EXAMPLE VALIDES
accession_num = folder.get_accession_number(DATA)
csv_file = f"{new_folder_example}/{accession_num}.csv"
seed = fastaReader.read_multi_fasta(DATA)
info_seq_dico = np.load(f"{DATA_SEQ_INFO}/{accession_num}.npy", allow_pickle="TRUE").item()
pid_file = np.load(f"{DATA_PID}/{accession_num}.pid.npy", allow_pickle="TRUE").item()
function_ex_save.ex_save(seed, info_seq_dico, pid_file, L,
                         csv_file, PID_INF, PID_SUP)


"""
DATA PRE-PROCESSING:
bash main_pid.sh > pid.out 2>&1 &
"""

# IMPORTS
import sys
import os
from pathlib import Path
import argparse
import numpy as np

import treatment.function_pid as function_pid
import utils.folder as folder
import utils.fastaReader as fastaReader

file = Path(__file__).resolve()
sys.path.append(file.parents[0])



# PARAMETERS

parser = argparse.ArgumentParser()
parser.add_argument("path_fasta_file",
                     help="path of file to describe the seq info from")
args = parser.parse_args()


DATA_PID = f"{file.parents[2]}/EXP_DIMANCHE/PID"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


# general folder managment
LIST_FOLDER = [DATA_PID]
for path in LIST_FOLDER:
    # if not os.path.exists(path):
    os.makedirs(path, exist_ok=True)

accession_num = folder.get_accession_number(args.path_fasta_file)
seed = fastaReader.read_multi_fasta(args.path_fasta_file)
path_file_pid = f"{DATA_PID}/{accession_num}.pid"
pid_couple = {}
n_seq = len(seed)
len_align = int(len(seed[0][1]))
for i in range(n_seq):
    pid_couple[seed[i][0]] = {}
        
for i in range(n_seq):
    for j in range(i, n_seq):
        current_pid = function_pid.pid(seed[i][1], seed[j][1],
                          ALPHABET, len_align)
        pid_couple[seed[i][0]][seed[j][0]] = current_pid
        pid_couple[seed[j][0]][seed[i][0]] = current_pid
        
np.save(path_file_pid, pid_couple)

  
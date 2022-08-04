"""
Preprocessing of data selection
"""

# IMPORTS
import os
import os.path
import argparse

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import example.function_seq_info as function_seq_info
import utils.folder as folder
import utils.fastaReader as fastaReader

# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("path_fasta_file",
                     help="path of file to describe the seq info from")
parser.add_argument("train_test", help="'TRAIN' or 'TEST")
args = parser.parse_args()

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

## main folders
DATA = args.path_fasta_file
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/2_PRE_EXAMPLE_{args.train_test}"




# FOLDER MANAGEMENT
new_folder_dico_seq_info = f"{DATA_RESULT}/SEQ_INFO"

list_folder = [DATA_RESULT, new_folder_dico_seq_info]
for path_folder in list_folder:
    os.makedirs(path_folder, exist_ok=True)



# SELECTION EXAMPLE VALIDES
accession_num = folder.get_accession_number(DATA)
seed = fastaReader.read_multi_fasta(DATA)
function_seq_info.seq_info(seed, accession_num, ALPHABET,
                               new_folder_dico_seq_info)

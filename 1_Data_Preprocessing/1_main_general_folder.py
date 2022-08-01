"""
FOLDER GENERAL PRE-PROCESSING
nohup python3 1_main_general_folder.py > 1.out 2>&1
"""

# IMPORTS
import sys
import os
from pathlib import Path
import time
from datetime import datetime

import treatment.stockholm as stockholm
import treatment.capitalizer as capitalizer
import treatment.description as description


file = Path(__file__).resolve()
sys.path.append(file.parents[0])



# PARAMETERS

# NAME_FOLDER_SOURCE = "MNHN_RESULT"
NAME_FOLDER_SOURCE = "MNHN_RESULT_FINAL"

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]



## main folders

FOLDER_SOURCE =  f"{file.parents[2]}/{NAME_FOLDER_SOURCE}"
DATA =  f"{FOLDER_SOURCE}/1_DATA"

NAME_MULTI_STOCKHOLM_FILE = "Pfam-A.seed"
NAME_MONO_STOCKHOLM_FOLDER = "Pfam_Stockholm"
NAME_FASTA_FOLDER = "Pfam_FASTA"
NAME_FASTA_FOLDER_UPPER = "Pfam_Upper"




# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

start = time.time()


## general folder managment
for path in [FOLDER_SOURCE, DATA]:
    os.makedirs(path, exist_ok='True')


print("\nMULTI --> MONO STOCKHOLM:")
path_file_multi_stockholm = f"{FOLDER_SOURCE}/{NAME_MULTI_STOCKHOLM_FILE}"
path_folder_mono_stockholm = f"{DATA}/{NAME_MONO_STOCKHOLM_FOLDER}"
stockholm.stockholm_separator(path_file_multi_stockholm, path_folder_mono_stockholm)


print("\nCONVERSION FROM STOCKHOLM FORMAT TO FASTA FORMAT:")
path_folder_fasta = f"{DATA}/{NAME_FASTA_FOLDER}"
stockholm.multi_stockholm_to_fasta(path_folder_mono_stockholm, path_folder_fasta)

path_character_percentage = f"{DATA}/character_stockholm.npy"
path_character_included_percentage = f"{DATA}/character_included_stockholm.npy"

description.data_count(path_folder_fasta, ALPHABET,
                       path_character_percentage,
                       path_character_included_percentage)

description.bar_plot_data_description(path_folder_fasta,
                                    path_character_percentage, "ALL CHARACTERS")
description.bar_plot_data_description(path_folder_fasta,
                    path_character_included_percentage , "STANDARD AMINO-ACIDS")



print("\nUPPER CASE:")
path_folder_fasta_upper = f"{DATA}/{NAME_FASTA_FOLDER_UPPER}"
capitalizer.multi_capitalization(path_folder_fasta, path_folder_fasta_upper)

end = time.time()
print(f"\nDONE IN: {'{:_}'.format(round(end - start, 4))} s")

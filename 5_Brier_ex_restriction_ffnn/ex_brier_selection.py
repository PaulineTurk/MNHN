# IMPORTS

import os
import os.path
import argparse
import time
import pandas as pd
import csv
import numpy as np
import random

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import utils.folder as folder

DATA_GENERAL = ""
DATA = ""
FRAC_SEED = ""

N_EX = 10_000_000



files = [x for x in Path(DATA).iterdir()]


with open(FRAC_SEED, newline='') as csvfrac:  # vérif que c bien read

    for file in files:
        accession_num = folder.get_accession_number(file)
        print(f"accession_num: {accession_num}")
        csvfrac[]

        

        with open(file, newline='') as csvfile:
            reader = csv.DictReader(csvfile)

            n_examples = len(reader)


            index_examples = random.choices()
            print(file)

            # counter_row = 0
            for row in reader:
                for direction in LIST_DIRECTON:
                    for position in LIST_ABS_POSITION:
                        dict_dict_3D[direction, position][row['aa_origin']][row['aa_destination']][row[f'aa_{direction}_{position}']] += 1
                # counter_row += 1
                # print(f"num_ex: {counter_row}")
        counter += 1
        print(round(100*counter/n_files, 2))
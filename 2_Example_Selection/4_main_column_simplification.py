
# IMPORTS

import os
import os.path
import argparse
import pandas as pd
import csv

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import utils.folder as folder


# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("path_csv_examples",
                     help="path of file with valid examples", type=str)
args = parser.parse_args()


L = 6
PID_INF = 40
PID_SUP = 50

DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/2_EXAMPLES"
new_folder_example = f"{DATA_RESULT}/EXAMPLES_{L}_{PID_INF}_{PID_SUP}_LIGHT"
os.makedirs(new_folder_example, exist_ok=True)

accession_num = folder.get_accession_number(args.path_csv_examples)


with open(args.path_csv_examples,"r") as source:
    rdr= csv.reader( source )
    with open(f"{new_folder_example}/{accession_num}.csv","w") as result:
        wtr= csv.writer( result )
        for r in rdr:
            # re-write all except the 4 first columns
            wtr.writerow(tuple(r[i] for i in range(4,len(r))))
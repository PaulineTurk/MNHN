"""SELECTION EX FOR BRIER EXPERIMENTS
nohup python3 ex_brier_selection.py > selection_ex_brier.out 2>&1 &
"""

# IMPORTS
import time
import csv
import numpy as np
import random
import math
import pandas as pd
from itertools import groupby

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import utils.folder as folder

DATA = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/EXAMPLES_6_40_50"
PATH_DICT_FRAC_SEED =  f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/frac_ex.csv"


start_dico = time.time()
# convert to a dico (do it once for all?)
dict_fraction = {}
with open(PATH_DICT_FRAC_SEED, newline='') as csv_frac:
    reader = csv.DictReader(csv_frac)
    for row in reader:
        dict_fraction[row["id_seed"]] = (row["num_ex"],row["frac_ex"])
end_dico = time.time()
print("DICO")
print(f"DONE IN: {round(end_dico - start_dico, 2)} s")


n_files = len(dict_fraction)

N_EX_TOTAL = 1000
L = 6

FILE_EXAMPLES = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/EX_BRIER_TRAIN.csv"

start = time.time()

files = [x for x in Path(DATA).iterdir()]


counter = 0

with open(FILE_EXAMPLES, 'w',
            encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header_context_ol = [f"aa_ol_{i}" for i in range(1, L+1)]
    header_context_or = [f"aa_or_{i}" for i in range(1, L+1)]
    header_context_dl = [f"aa_dl_{i}" for i in range(1, L+1)]
    header_context_dr = [f"aa_dr_{i}" for i in range(1, L+1)]
    header_info = ["counter", "pid", "name_origin", "name_destination", "aa_origin", "aa_destination"]

    header = header_info + header_context_ol + header_context_or + header_context_dl + header_context_dr
    writer.writerow(header)


    for file in files:
        accession_num = folder.get_accession_number(file)
        n_ex_in_seed_total = int(dict_fraction[accession_num][0])
        if n_ex_in_seed_total != 0:
            fraction = float(dict_fraction[accession_num][1])
            N_EX_SEED_FLOAT = fraction*N_EX_TOTAL
            decimal_part, integer_part = math.modf(float(N_EX_SEED_FLOAT))
            proba = random.uniform(0, 1)
            N_EX_SEED_INT = int(integer_part + int(proba < decimal_part))


            with open(file, newline='') as csvfile:
                reader = csvfile.read().splitlines()

                # n_lines pris avant
                list_index = random.choices([i for i in range(1, n_ex_in_seed_total+1)], k=N_EX_SEED_INT) # first line headers
                list_index.sort()
                list_index_organised = [[k,len(list(v))] for k,v in groupby(list_index)]
                if list_index_organised != []:
        
                    i = 0 # index in final_list_index
                    count = 0 # count of examples selected in the current seed
                    index_row = 0
                    while count < N_EX_SEED_INT: # strictement ?
                        # row = next(csvfile)
                        index_row += 1
                        if index_row == list_index_organised[i][0]:
                            data = reader[index_row].split(",")
                            # data = str(row).split(",")
                            # print(data)
                            for j in range(list_index_organised[i][1]):
                                writer.writerow(data)
                            count += list_index_organised[i][1]
                            i += 1
                    

        counter += 1
        print(round(100*counter/n_files, 2))

end = time.time()

print("SELECTION EXAMPLES")
print(f"DONE IN: {round(end - start, 2)} s")

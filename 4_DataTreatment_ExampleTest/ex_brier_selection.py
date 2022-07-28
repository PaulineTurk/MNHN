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

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import utils.folder as folder

DATA = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/test_concat"
PATH_DICT_FRAC_SEED =  f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/frac_ex.csv"


start_dico = time.time()
# convert to a dico (do it once for all?)
dict_fraction = {}
with open(PATH_DICT_FRAC_SEED, newline='') as csv_frac:
    reader = csv.DictReader(csv_frac)
    for row in reader:
        dict_fraction[row["id_seed"]] = row["frac_ex"]
end_dico = time.time()
print("DICO")
print(f"DONE IN: {round(end_dico - start_dico, 2)} s")


n_files = len(dict_fraction)

N_EX_TOTAL = 1_000
L = 6

FILE_EXAMPLES = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES/EX_BRIER_TRAIN_mini_fast_old.csv"

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
    header_info = ["pid", "name_origin", "name_destination", "aa_origin", "aa_destination"]

    header = header_info + header_context_ol + header_context_or + header_context_dl + header_context_dr
    writer.writerow(header)


    for file in files:
        accession_num = folder.get_accession_number(file)
        print(f"accession_num: {accession_num}")
        fraction = float(dict_fraction[accession_num])
        
        if fraction != 0:
            N_EX_SEED_FLOAT = fraction*N_EX_TOTAL
            decimal_part, integer_part = math.modf(float(N_EX_SEED_FLOAT))
            proba = random.uniform(0, 1)
            N_EX_SEED_INT = int(integer_part + int(proba < decimal_part))

            with open(file, newline='') as csvfile:
                reader = csv.reader(csvfile)
                lines = [tuple(line) for line in reader]
                # print(f"lines:  {lines}")
                chosen_rows = random.choices(lines, k=N_EX_SEED_INT)   # header peut etre selectionné, c pas bon ca !
            #     # print(f"chosen rows: {chosen_rows}")
                np.random.choice([1, 2, 3], size=10, replace=True)
                for data in chosen_rows:
                    writer.writerow(data)


            # with open(file, newline='') as csvfile:
            #     df = pd.read_csv(csvfile)
            #     n_lines = len(df)  # voir avec header
            #     list_index = [i for i in range(n_lines-1)]  # verif headers ...
            #     index_selected = random.choices(list_index, k = N_EX_SEED_INT)
            #     for i in index_selected:
            #         writer.writerow(tuple(df.iloc[i]))

        counter += 1
        print(round(100*counter/n_files, 2))

end = time.time()

print("SELECTION EXAMPLES")
print(f"DONE IN: {round(end - start, 2)} s")

"""
2D_VERSIONS
nohup python3 2_main_2d_versions.py > 2.out 2>&1 &
"""

# IMPORTS
import numpy as np
import time
from datetime import datetime

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import table2d.function_table2d as function_table2d


# PARAMETERS

PID_INF = 40
PID_SUP = 50
L = 6

DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_2D/{L}_{PID_INF}_{PID_SUP}"

LIST_PSEUDO_COUNTER_2D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10,  0),
                          pow(10,  1),
                          pow(10,  2),
                          pow(10,  3)]

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]






# PROGRAM
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)

start_total = time.time()
start = time.time()

print("")
print(f"count_1d")
print("")
path_count_2d = f"{DATA_RESULT}/count_2d.npy"
dict_count_2d = np.load(path_count_2d, allow_pickle='True').item()

# initialisation dict_count_1d
dict_count_1d = {}
for aa in ALPHABET:
    dict_count_1d[aa] = 0

total_aa = 0
for aa_origin in ALPHABET:
    for aa_destination in ALPHABET:
        dict_count_1d[aa_origin] += dict_count_2d[aa_origin][aa_destination]
        total_aa += dict_count_2d[aa_origin][aa_destination]
end = time.time()
print(f"\nDONE IN: {round(end - start, 2)} s")
np.save(f"{DATA_RESULT}/count_1d.npy", dict_count_1d)  # verif total
                                                       # verif nom .npy.npy


start = time.time()
print("")
print(f"freq_1d")
print("")
dict_freq_1d = {}
for aa in ALPHABET:
    dict_freq_1d[aa] = dict_count_1d[aa]/total_aa
end = time.time()
print(f"\nDONE IN: {round(end - start, 2)} s")
np.save(f"{DATA_RESULT}/freq_1d.npy", dict_freq_1d)


start = time.time()
print("")
print(f"freq_2d")
print("")
dict_freq_2d = {}
for aa_origin in ALPHABET:
    dict_freq_2d[aa_origin] = {}
    for aa_destination in ALPHABET:
        dict_freq_2d[aa_origin][aa_destination] = dict_count_2d[aa_origin][aa_destination]/total_aa
end = time.time()
print(f"\nDONE IN: {round(end - start, 2)} s")
np.save(f"{DATA_RESULT}/freq_2d.npy", dict_freq_2d)



print("")
print(f"proba_pseudo_count")
print("")
for pseudo_counter_2d in LIST_PSEUDO_COUNTER_2D:
    pseudo_counter_2d = float(pseudo_counter_2d)
    print("")
    print(pseudo_counter_2d)
    print("")
    path_freq_AA = f"{DATA_RESULT}/freq_1d.npy"           # dict_freq_1d
    path_freq_AA_couple = f"{DATA_RESULT}/freq_2d.npy"    # dict_freq_2d
    path_file_Result = f"{DATA_RESULT}/proba_{pseudo_counter_2d}.npy" # dict_proba_pseudo_counter_2d                  

    function_table2d.proba_conditional_weighted(path_freq_AA, path_freq_AA_couple,
                                                pseudo_counter_2d,
                                                ALPHABET,
                                                path_file_Result)

print("")
print(f"score")
print("")
path_file_Result = f"{DATA_RESULT}/score.npy"  # dict_score
function_table2d.score(path_freq_AA,
                      path_freq_AA_couple,
                      path_file_Result,
                      ALPHABET,
                      scale_factor=2)

end = time.time()
print(f"\nDONE IN: {round(end - start_total, 2)} s")

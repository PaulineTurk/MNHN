
# IMPORTS
import numpy as np
import time

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])



# FUNCTIONS
def table_3d_proba(table_3d_count, alphabet):

    intra_couple_count = {}
    for aa_o in alphabet:
        intra_couple_count[aa_o] = {}
        for aa_d in alphabet:
            intra_couple_count[aa_o][aa_d] = 0
            for aa_c in alphabet:
                intra_couple_count[aa_o][aa_d] += table_3d_count[aa_o][aa_d][aa_c]

    table_3d_proba = {}
    for aa_o in alphabet:
        table_3d_proba[aa_o] = {}
        for aa_d in alphabet:
            table_3d_proba[aa_o][aa_d] = {}
            for aa_c in alphabet:
                count_triplet = table_3d_count[aa_o][aa_d][aa_c]
                value_intra_couple = intra_couple_count[aa_o][aa_d]
                if value_intra_couple != 0:
                    table_3d_proba[aa_o][aa_d][aa_c] = count_triplet/value_intra_couple
                else:
                    table_3d_proba[aa_o][aa_d][aa_c] = 0
    return table_3d_proba


def normalisation(table_3d_proba, alphabet):
    table_3d_proba_normal = {}
    for aa_o in alphabet:
        table_3d_proba_normal[aa_o] = {}
        for aa_d in alphabet:
            table_3d_proba_normal[aa_o][aa_d] = {}
            sum_line = sum([table_3d_proba[aa_o][aa_d][aa_c] for aa_c in alphabet])
            for aa_c in alphabet:
                if sum_line != 0:
                    table_3d_proba_normal[aa_o][aa_d][aa_c] = table_3d_proba[aa_o][aa_d][aa_c]/sum_line
    return table_3d_proba_normal


def proba_normal_weighting(path_folder_table_3d_count,
                           path_folder_table_2d_proba,
                           pid_inf, pid_sup,
                           pseudo_counter_2d, pseudo_counter_3d,
                           relative_index_context, RefSeq,
                           weight, alphabet,
                           path_folder_save):
    # STR CONVERSION
    pid_inf = str(pid_inf)
    pid_sup = str(pid_sup)
    pseudo_counter_2d = str(pseudo_counter_2d)
    pseudo_counter_3d = str(pseudo_counter_3d)
    relative_index_context = str(relative_index_context)

    # 2D PROBA LOAD
    path_2d_proba = f"{path_folder_table_2d_proba}/{pid_inf}_{pid_sup}_{pseudo_counter_2d}/proba.npy"
    table_2d_proba = np.load(path_2d_proba, allow_pickle='TRUE').item()
    print(f"2D PROBA PATH: {path_2d_proba}")

    # 3D COUNT LOAD
    path_3d_count = f"{path_folder_table_3d_count}/{pid_inf}_{pid_sup}_{pseudo_counter_3d}/{relative_index_context}_{RefSeq}.npy"
    table_3d_count = np.load(path_3d_count, allow_pickle='TRUE').item()
    print(f"3D COUNT PATH: {path_3d_count}")

    # 3D COUNT TO PROBA (not weighted)
    table_3d_proba_init = table_3d_proba(table_3d_count, alphabet)

    # 3D PROBA WEIGHTING
    start = time.time()
    table_3d_proba_w = {}
    for aa_o in alphabet:
        table_3d_proba_w[aa_o] = {}
        for aa_d in alphabet:
            table_3d_proba_w[aa_o][aa_d] = {}
            for aa_c in alphabet:
                context_estimation = table_3d_proba_init[aa_o][aa_d][aa_c]
                context_free = table_2d_proba[aa_o][aa_d]
                table_3d_proba_w[aa_o][aa_d][aa_c] = (context_estimation + weight*context_free)/ (1+weight)
    end = time.time()
    diff = end - start
    print(f"TIME WEIGHTING : {diff:.2f}s")

    # 3D PROBA NORMALISATION
    start = time.time()
    table_3d_proba_w_n = normalisation(table_3d_proba_w, alphabet)
    end = time.time()
    diff = end - start
    print(f"TIME NORMALISATION : {diff:.2f}s")

    # SAVING
    path_file_save = f"{path_folder_save}/{pid_inf}_{pid_sup}_{pseudo_counter_3d}_{str(weight)}/{relative_index_context}_{RefSeq}"
    np.save(path_file_save, table_3d_proba_w_n)
    print(f"PATH 3D_PROBA WEIGHTED AND NORMALISED: {path_file_save}.npy")

    return table_3d_proba_w_n



def sum_line(table_3d_proba, alphabet, aa_o, aa_d):
    """
    Sum of the conditional probabilities on each line (aa_o and aa_d fixed)
    must be equal to 1 to respect the total probability formula.
    """
    sum_line = 0
    for aa_c in alphabet:
        sum_line += table_3d_proba[aa_o][aa_d][aa_c]
    print(sum_line)
    return sum_line



def sum_plate(table_3d_proba):
    """
    cond_proba: cube of the conditional probabilities of each valid triplet.

    The sum of the conditional probabilities on each horizontal level of the cube
    must be equal to the length of an edge of the cube to respect the total probability formula.
    """
    for aa_o in table_3d_proba:
        sum_plate = 0
        for aa_d in table_3d_proba[aa_o]:
            for aa_c in table_3d_proba[aa_o][aa_d]:
                sum_plate += table_3d_proba[aa_o][aa_d][aa_c]
        print(f"{aa_o}, {sum_plate}")

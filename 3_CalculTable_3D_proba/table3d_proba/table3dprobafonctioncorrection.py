
# IMPORTS
import numpy as np
import time

import sys  
from pathlib import Path

from pkg_resources import VersionConflict 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])



# FUNCTIONS


# à tester
def freq_3d(path_table_3d_count,
            path_save_table_3d_freq,
            alphabet):
    print(f"path table_3d_count: {path_table_3d_count}")
    table_3d_count = np.load(path_table_3d_count, allow_pickle='TRUE').item()
    # N_EXAMPLE_TRAIN COUNT
    n_example_train = 0
    for aa_origine in alphabet:
        for aa_destination in alphabet:
            for aa_context in alphabet:
                n_example_train += table_3d_count[aa_origine][aa_destination][aa_context]
    print(f"N_EX_TRAIN: {'{:_}'.format(n_example_train)}")

    # FROM table_3d_count to table_3d_freq
    table_3d_freq = {}
    for aa_origine in alphabet:
        table_3d_freq[aa_origine] = {}
        for aa_destination in alphabet:
            table_3d_freq[aa_origine][aa_destination] = {}
            for aa_context in alphabet:
                table_3d_freq[aa_origine][aa_destination][aa_context] = table_3d_count[aa_origine][aa_destination][aa_context]/n_example_train

    np.save(path_save_table_3d_freq, table_3d_freq)

def sum_freq_3D(path_table_3d_freq, alphabet):
    table_3d_freq = np.load(path_table_3d_freq, allow_pickle='TRUE').item()
    sum_freq = 0
    for aa_origine in alphabet:
        for aa_destination in alphabet:
            for aa_context in alphabet:
                sum_freq += table_3d_freq[aa_origine][aa_destination][aa_context]
    print(f"SUM_FREQ: {sum_freq}")



# à tester
def proba_3d(path_table_3d_freq,
             path_table_2d_proba_pc,
             path_table_1d_freq,
             pseudo_counter_3d,
             alphabet,
             path_save_table_3d_pc):

    table_3d_freq = np.load(path_table_3d_freq, allow_pickle='TRUE').item()
    table_1d_freq = np.load(path_table_1d_freq, allow_pickle='TRUE').item()
    table_2d_proba_pc = np.load(path_table_2d_proba_pc, allow_pickle='TRUE').item()

    table_proba_3d_pc = {}
    for aa_origine in alphabet:
        table_proba_3d_pc[aa_origine] = {}
        for aa_destination in alphabet:
            table_proba_3d_pc[aa_origine][aa_destination] = {}
            vector_context = []
            proba_2d_pc = table_2d_proba_pc[aa_origine][aa_destination]
            for aa_context in alphabet:
                freq_triplet = table_3d_freq[aa_origine][aa_destination][aa_context]
                freq_context = table_1d_freq[aa_context]
                freq_origine = table_1d_freq[aa_origine]
                vector_context.append(freq_triplet + pseudo_counter_3d*proba_2d_pc*freq_origine*freq_context)
            denominator = sum(vector_context)

            for index, aa_context in enumerate(alphabet):
                if denominator != 0:
                    table_proba_3d_pc[aa_origine][aa_destination][aa_context] = vector_context[index]/denominator
                else:
                    table_proba_3d_pc[aa_origine][aa_destination][aa_context] = 0
                    

    np.save(path_save_table_3d_pc, table_proba_3d_pc)







# def table_3d_proba(table_3d_count, alphabet):

#     intra_couple_count = {}
#     for aa_o in alphabet:
#         intra_couple_count[aa_o] = {}
#         for aa_d in alphabet:
#             intra_couple_count[aa_o][aa_d] = 0
#             for aa_c in alphabet:
#                 intra_couple_count[aa_o][aa_d] += table_3d_count[aa_o][aa_d][aa_c]

#     table_3d_proba = {}
#     for aa_o in alphabet:
#         table_3d_proba[aa_o] = {}
#         for aa_d in alphabet:
#             table_3d_proba[aa_o][aa_d] = {}
#             for aa_c in alphabet:
#                 count_triplet = table_3d_count[aa_o][aa_d][aa_c]
#                 value_intra_couple = intra_couple_count[aa_o][aa_d]
#                 if value_intra_couple != 0:
#                     table_3d_proba[aa_o][aa_d][aa_c] = count_triplet/value_intra_couple
#                 else:
#                     table_3d_proba[aa_o][aa_d][aa_c] = 0
#     return table_3d_proba


# def normalisation(table_3d_proba, alphabet):
#     table_3d_proba_normal = {}
#     for aa_o in alphabet:
#         table_3d_proba_normal[aa_o] = {}
#         for aa_d in alphabet:
#             table_3d_proba_normal[aa_o][aa_d] = {}
#             sum_line = sum([table_3d_proba[aa_o][aa_d][aa_c] for aa_c in alphabet])
#             for aa_c in alphabet:
#                 if sum_line != 0:
#                     table_3d_proba_normal[aa_o][aa_d][aa_c] = table_3d_proba[aa_o][aa_d][aa_c]/sum_line
#     return table_3d_proba_normal


# def proba_normal_weighting(path_folder_table_3d_count,
#                            path_folder_table_2d_proba,
#                            pid_inf, pid_sup,
#                            pseudo_counter_2d, pseudo_counter_3d,
#                            relative_index_context, RefSeq,
#                            weight, alphabet,
#                            path_folder_save):
#     # STR CONVERSION
#     pid_inf = str(pid_inf)
#     pid_sup = str(pid_sup)
#     pseudo_counter_2d = str(pseudo_counter_2d)
#     pseudo_counter_3d = str(pseudo_counter_3d)
#     relative_index_context = str(relative_index_context)

#     # 2D PROBA LOAD
#     path_2d_proba = f"{path_folder_table_2d_proba}/{pid_inf}_{pid_sup}_{pseudo_counter_2d}/proba.npy"
#     table_2d_proba = np.load(path_2d_proba, allow_pickle='TRUE').item()
#     print(f"2D PROBA PATH: {path_2d_proba}")

#     # 3D COUNT LOAD
#     path_3d_count = f"{path_folder_table_3d_count}/{pid_inf}_{pid_sup}_{pseudo_counter_3d}/{relative_index_context}_{RefSeq}.npy"
#     table_3d_count = np.load(path_3d_count, allow_pickle='TRUE').item()
#     print(f"3D COUNT PATH: {path_3d_count}")

#     # 3D COUNT TO PROBA (not weighted)
#     table_3d_proba_init = table_3d_proba(table_3d_count, alphabet)

#     # 3D PROBA WEIGHTING
#     start = time.time()
#     table_3d_proba_w = {}
#     for aa_o in alphabet:
#         table_3d_proba_w[aa_o] = {}
#         for aa_d in alphabet:
#             table_3d_proba_w[aa_o][aa_d] = {}
#             for aa_c in alphabet:
#                 context_estimation = table_3d_proba_init[aa_o][aa_d][aa_c]
#                 context_free = table_2d_proba[aa_o][aa_d]
#                 table_3d_proba_w[aa_o][aa_d][aa_c] = (context_estimation + weight*context_free)/ (1+weight)
#     end = time.time()
#     diff = end - start
#     print(f"TIME WEIGHTING : {diff:.2f}s")

#     # 3D PROBA NORMALISATION
#     start = time.time()
#     table_3d_proba_w_n = normalisation(table_3d_proba_w, alphabet)
#     end = time.time()
#     diff = end - start
#     print(f"TIME NORMALISATION : {diff:.2f}s")

#     # SAVING
#     path_file_save = f"{path_folder_save}/{pid_inf}_{pid_sup}_{pseudo_counter_3d}_{str(weight)}/{relative_index_context}_{RefSeq}"
#     np.save(path_file_save, table_3d_proba_w_n)
#     print(f"PATH 3D_PROBA WEIGHTED AND NORMALISED: {path_file_save}.npy")

#     return table_3d_proba_w_n



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

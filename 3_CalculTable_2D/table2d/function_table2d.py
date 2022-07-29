"""
FUNCTIONS FOR NON-CONTEXTUAL INFORMATION
"""

# IMPORTS

import pandas as pd
from math import log2
import numpy as np
import blosum as bl
import seaborn as sb
import matplotlib.pyplot as plt
import time


import sys  
from pathlib import Path  
file = Path(__file__).resolve()
sys.path.append(file.parents[1])




def proba_conditional_weighted(path_freq_AA: str, path_freq_AA_couple: str,
                               pseudo_counter_2d: float,
                               alphabet: list,
                               path_file_Result: str):
    start = time.time()
    # LOAD
    dico_freq_AA = np.load(path_freq_AA, allow_pickle=True).item()
    dico_freq_AA_couple = np.load(path_freq_AA_couple, allow_pickle=True).item()

    dico_conditional_proba_weight = {}
    for char_1 in alphabet:
        dico_conditional_proba_weight[char_1] = {}
        for char_2 in alphabet:
            numerator = (dico_freq_AA_couple[char_1][char_2]
                        + pseudo_counter_2d*dico_freq_AA[char_1]*dico_freq_AA[char_2])
            denumerator = (1 + pseudo_counter_2d)*dico_freq_AA[char_1]
            if denumerator != 0:
                dico_conditional_proba_weight[char_1][char_2] = numerator/denumerator
            else:
                dico_conditional_proba_weight[char_1][char_2] = 0

    np.save(path_file_Result, dico_conditional_proba_weight)
    end = time.time()
    print(f"2D PROBA: time {'{:_}'.format(round(end - start, 4))}s")




def score(path_freq_AA, path_freq_AA_couple,
          path_file_Result, alphabet, scale_factor=2):
    """
    Compute and save table_2d_score
    """
    start = time.time()

    # LOAD: FREQ DICTIONARIES
    dico_freq_AA = np.load(path_freq_AA, allow_pickle=True).item()
    dico_freq_AA_couple = np.load(path_freq_AA_couple, allow_pickle=True).item()

    table_2d = {}
    for aa_1 in alphabet:
        table_2d[aa_1] = {}
        for aa_2 in alphabet:
            if dico_freq_AA_couple[aa_1][aa_2] != 0:
                table_2d[aa_1][aa_2] = round(scale_factor * log2(dico_freq_AA_couple[aa_1][aa_2]
                                            / (dico_freq_AA[aa_1] * dico_freq_AA[aa_2])))
            else:
                table_2d[aa_1][aa_2] = 0
    end = time.time()
    print(f"2D SCORE: time {'{:_}'.format(round(end - start, 4))} s")

    np.save(path_file_Result, table_2d)




def table_2d_heatmap(matrix, path_folder_Result, title, size_annot = 3):
    """
    Save the heatmap of the matrix in path_folder_Result
    """
    heatmap_matrix = np.transpose(pd.DataFrame.from_dict(matrix))
    heatmap = sb.heatmap(heatmap_matrix, annot = True, annot_kws = {"size": size_annot}, fmt = '.2g')
    plt.yticks(rotation=0)
    heatmap_figure = heatmap.get_figure()
    plt.title(title, loc='center', wrap=True)
    plt.close()
    path_save_fig = f"{path_folder_Result}/{title}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)



def table_2d_visualisation(table_2d):
    """
    Visualisation of the matrix
    """
    df_table_2d = np.transpose(pd.DataFrame.from_dict(table_2d))
    print(df_table_2d)


def sum_line(table_2d):
    """
    To check that the sum of a line is equal to one
    for the conditional probability matrix
    """
    df_table_2d= np.transpose(pd.DataFrame.from_dict(table_2d))
    sum_line = df_table_2d.sum(axis=1)
    print("SUM OVER A CHARACTER OF DESTINATION:")
    print("")
    print(sum_line)


def table_2d_difference(table_2d, alphabet, pid_inf_ref):
    """
    Quantify the distance between the table_2d_score computed and a blosum of reference
    """
    # IMPORTS
    blosum_ref = bl.BLOSUM(pid_inf_ref)

    # INITIALISATION
    matrix_diff = {}
    difference = 0
    count = 0

    # DIFFERENCES EVALUATION
    for aa1 in alphabet:
        matrix_diff[aa1] = {}
        for aa2 in alphabet:
            matrix_diff[aa1][aa2] = int(table_2d[aa1][aa2] - blosum_ref[aa1 + aa2])
            difference += matrix_diff[aa1][aa2]
            count += 1
    average_difference  = round(difference/count, 2)

    return matrix_diff, pid_inf_ref, average_difference

def min_2D(dic_2D, alphabet):
    minimum = 1
    for char in alphabet:
        minimum_temp = min(dic_2D[char].values())
        print(f"min_1: {minimum_temp}")
        if float(minimum_temp) < minimum:
            minimum = minimum_temp
    return minimum

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
from tqdm import tqdm
import time


import sys  
from pathlib import Path  
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.fastaReader as fastaReader
import utils.folder as folder


# FUNCTIONS
def count_for_table_2d(num_accession: str, path_folder_pid: str,
                       seed, pid_inf: int, pid_sup: int,
                       dico_count_AA: dict,
                       dico_count_AA_couple: dict,
                       alphabet: list):
    """Count the number of each character of the alphabet and the each
       couple of characters in one multi-sequence alignment (seed.
       The counts are concatenated to the previous counts on other seeds.

    Args:
        num_accession (str): accession to a multi-sequence alignment
        path_folder_pid (str): path to the percentage of identity information
        seed: a multi-sequence alignment (name, seq)
        pid_inf (float): a pairwise alignment must have a pid >= pid_inf
        pid_sup (float): a pairwise alignment must have a pid < pid_inf
        dico_count_AA (dict): count of each character in alphabet
        dico_count_AA_couple (dict): count of each couple of character in alphabet
        alphabet (list): alphabet (list): list of character included

    Returns:
        dict, dict: dico_count_AA, dico_count_AA_couple
    """

    pid_couple = np.load(f"{path_folder_pid}/{num_accession}.pid.npy", allow_pickle='TRUE').item()
    nb_seq = len(seed)
    for i in range(nb_seq -1):
        name_1, seq_1 = seed[i]
        for j in range(i + 1, nb_seq):
            name_2 ,seq_2 = seed[j]
            if pid_couple[name_1][name_2] >= pid_inf and pid_couple[name_1][name_2] < pid_sup:
                for (aa_1, aa_2) in zip(seq_1, seq_2):
                    if aa_1 in alphabet and aa_2 in alphabet:

                        # count AA
                        dico_count_AA[aa_1] += 2
                        dico_count_AA[aa_2] += 2

                        # count couple AA
                        if aa_1 == aa_2:
                            dico_count_AA_couple[aa_1][aa_2] += 2
                        else:
                            dico_count_AA_couple[aa_1][aa_2] += 1
                            dico_count_AA_couple[aa_2][aa_1] += 1

    return dico_count_AA, dico_count_AA_couple



def multi_count_for_table_2d(path_folder_fasta: str, path_folder_pid: str,
                           alphabet: list, pid_inf: int, pid_sup: int,
                           path_folder_Result: str):
    """
    Iterate count_for_table_2d on the files included in path_folder_fasta.
    """

    # INITIALISATION: AMINO-ACID COUPLE COUNT
    dico_count_AA_couple = {}
    for aa_1 in alphabet:
        dico_count_AA_couple[aa_1] = {}
        for aa_2 in alphabet:
            dico_count_AA_couple[aa_1][aa_2] = 0

    # INITIALISATION: AMINO-ACID COUNT
    dico_count_AA = {}
    for aa in alphabet:
        dico_count_AA[aa] = 0

    files = [x for x in Path(path_folder_fasta).iterdir()]
    n_files = len(files)

    # COUNT
    start = time.time()
    for file_counter in tqdm(range(n_files), desc = "2d count",
                             ncols= 100, mininterval=60):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        seed_train = fastaReader.read_multi_fasta(file)
        dico_count_AA, dico_count_AA_couple = count_for_table_2d(
                                                accession_num,
                                                path_folder_pid,
                                                seed_train,
                                                pid_inf, pid_sup,
                                                dico_count_AA,
                                                dico_count_AA_couple,
                                                alphabet)

    end = time.time()
    diff = end - start
    items_per_second = n_files/diff
    print(f"2D COUNT: {'{:_}'.format(round(end - start, 4))}s | {items_per_second:.2f}it/s")

    # SAVE
    path_count_AA = f"{path_folder_Result}/count_1d"
    np.save(path_count_AA, dico_count_AA)

    path_count_AA_couple = f"{path_folder_Result}/count_2d"
    np.save(path_count_AA_couple, dico_count_AA_couple)

    # PARAMETERS NOT EVALUATED CHECK
    n_parameters_not_evaluated = 0
    for aa_1 in alphabet:
        for aa_2 in alphabet:
            if dico_count_AA_couple[aa_1][aa_2] == 0:
                n_parameters_not_evaluated += 1
    percentage_not_evaluated = round(100*n_parameters_not_evaluated/len(alphabet)**2)
    print(f"\nPARAMETERS_NOT_EVALUATED: {percentage_not_evaluated} %")

    return percentage_not_evaluated




def proba_conditional_weighted(path_freq_AA: str, path_freq_AA_couple: str,
                               pseudo_counter_2d: float,
                               alphabet: list,
                               path_folder_Result: str):
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

    path_dico_conditional_proba_weight = f"{path_folder_Result}/proba_{pseudo_counter_2d}"
    np.save(path_dico_conditional_proba_weight, dico_conditional_proba_weight)
    end = time.time()
    print(f"2D PROBA: time {'{:_}'.format(round(end - start, 4))}s")


def freq(path_count_AA: str, path_count_AA_couple: str,
         alphabet: list, path_folder_Result: str):

    start = time.time()

    # LOAD: COUNT DICTIONARIES
    dico_count_AA = np.load(path_count_AA, allow_pickle=True).item()
    dico_count_AA_couple = np.load(path_count_AA_couple, allow_pickle=True).item()

    # count_AA to freq_AA
    denominator = sum(dico_count_AA.values())
    dico_freq_AA = {key: value/denominator for key, value in dico_count_AA.items()}

    # SAVE
    path_dico_freq_AA = f"{path_folder_Result}/freq_1d"
    np.save(path_dico_freq_AA, dico_freq_AA)

    # count_AA_couple to freq_AA_couple
    n_ex_train = 0
    for aa_1 in alphabet:
        for aa_2 in alphabet:
            n_ex_train += dico_count_AA_couple[aa_1][aa_2]
    print(f"N_EXAMPLES_TRAIN: {'{:_}'.format(int(n_ex_train))}")

    dico_freq_AA_couple = {}
    for aa_1 in alphabet:
        dico_freq_AA_couple[aa_1] = {}
        for aa_2 in alphabet:
            dico_freq_AA_couple[aa_1][aa_2] = dico_count_AA_couple[aa_1][aa_2]/n_ex_train

    # SAVE
    path_dico_freq_AA_couple = f"{path_folder_Result}/freq_2d"
    np.save(path_dico_freq_AA_couple, dico_freq_AA_couple)

    end = time.time()
    print(f"2D FREQ: time {'{:_}'.format(round(end - start, 4))}s")



def score(path_freq_AA, path_freq_AA_couple,
        path_folder_Result, alphabet, scale_factor=2):
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

    path_matrix = f"{path_folder_Result}/score"
    np.save(path_matrix, table_2d)




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

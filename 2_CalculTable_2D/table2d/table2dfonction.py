"""
FUNCTIONS FOR NON-CONTEXTUAL INFORMATION
"""

# IMPORTS

import pandas as pd
from math import log2
import numpy as np
import os 
import blosum as bl
import seaborn as sb
import matplotlib.pyplot as plt
from tqdm import tqdm
import time


import sys  
from pathlib import Path  
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
from utils.timer import Timer
import utils.fastaReader as fastaReader
import utils.folder as folder


# FUNCTIONS
def count_for_table_2d(num_accession, path_folder_pid, seed, pid_inf, pid_sup,
                     count_AA, nb_AA, count_coupleAA, nb_coupleAA, nb_ex_train, list_residu):
    """
    In the seed with id num_accession, count the number of each valid amino acid
    and the number of valid couple
    
    num_accession: pid of the seed
    name_folder_pid: the pid file of this seed
    seed: the (name, seq) tuples of this seed
    pid_inf: smaller pid to validate a couple of sequence
    count_AA: dictionary of the count of each valid amino acid to cumulate this count on many seeds
    nb_AA: count of all valid amino acids
    count_coupleAA: dictionary of the count of each valid couple of amino acids to cumulate this count on many seeds
    nb_coupleAA: count of all valid couple of amino acids
    list_residu: list of valid amino acids
    """

    pid_couple = np.load(f"{path_folder_pid}/{num_accession}.pid.npy", allow_pickle='TRUE').item()
    nb_seq = len(seed)
    for i in range(nb_seq -1): # dans test préliminaire il n'y avait pas le -1
        name_1, seq_1 = seed[i]
        for j in range(i + 1, nb_seq):
            name_2 ,seq_2 = seed[j] 
            if pid_couple[name_1][name_2] >= pid_inf and pid_couple[name_1][name_2] < pid_sup:
                for (aa_1, aa_2) in zip(seq_1, seq_2):  
                    if aa_1 in list_residu and aa_2 in list_residu:
                        nb_ex_train += 2   # exemple selectionné et son symétrique

                        # count AA
                        count_AA[aa_1] += 1
                        count_AA[aa_2] += 1
                        nb_AA += 2

                        # count couple AA
                        if aa_1 == aa_2:
                            count_coupleAA[aa_1][aa_2] += 2
                        else:
                            count_coupleAA[aa_1][aa_2] += 1
                            count_coupleAA[aa_2][aa_1] += 1
                        nb_coupleAA += 2

    return count_AA, nb_AA, count_coupleAA, nb_coupleAA, nb_ex_train



def multi_count_for_table_2d(path_folder_fasta, path_folder_pid,
                           list_residu, pid_inf, pid_sup,
                           path_folder_Result,
                           pseudo_compte):
    """
    Iterate count_for_table_2d on the files included in path_folder_fasta.
    """

    nb_ex_train = 0
    nb_AA = 0
    nb_coupleAA = 0
    len_alphabet = len(list_residu)
 
    # intialisation avec pseudo-compte de la table_2d de comptage
    count_coupleAA = {}
    for aa_1 in list_residu:
        count_coupleAA[aa_1] = {}
        for aa_2 in list_residu:
            count_coupleAA[aa_1][aa_2] = pseudo_compte
            nb_coupleAA += pseudo_compte

    # intialisation avec pseudo-compte du compte de chaque acide aminé dans la table_2d de comptage
    count_AA = {}
    for aa in list_residu:
        count_AA[aa] = len_alphabet*pseudo_compte
        nb_AA += len_alphabet*pseudo_compte



    path, dirs, files = next(os.walk(path_folder_fasta))
    nb_files = len(files)

    # liste des PosixPath des alignements d'apprentissage
    files = [x for x in Path(path_folder_fasta).iterdir()]
    start = time.time()
    for file_counter in tqdm(range(nb_files), desc = "calcul de table_2d_count",
                             ncols= 100, mininterval=60):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        seed_train = fastaReader.read_multi_fasta(file)
        count_AA, nb_AA, count_coupleAA, nb_coupleAA, nb_ex_train = count_for_table_2d(accession_num, path_folder_pid, seed_train, pid_inf, pid_sup,
                                                                                count_AA, nb_AA, count_coupleAA, nb_coupleAA, nb_ex_train, list_residu)

    end = time.time()
    diff = end - start
    items_per_second = nb_files/diff
    print(f"Count of amino acid and couple of amino acids: {diff:.2f}s | {items_per_second:.2f}it/s")
    # print("\ncount_AA:\n", count_AA)
    # print("\nnb_AA:\n", nb_AA)
    # print("\ncount_coupleAA:\n", count_coupleAA)
    # print("\nnb_coupleAA:\n", nb_coupleAA)

    # vérification que chacun des 400 paramètres ont été évalué par au moins un exemple d'apprentissage
    count_couple_not_evaluated = 0
    for aa_1 in list_residu:
        for aa_2 in list_residu:
            if count_coupleAA[aa_1][aa_2] == pseudo_compte:
                count_couple_not_evaluated += 1
    percentage_couple_not_evaluated = 100*count_couple_not_evaluated/(len(list_residu)**2)
    print(f"Pourcentage de couples d'acides aminés non estimés: {percentage_couple_not_evaluated} %")

    print("Nombre d'exemples d'apprentissage :",  '{:_}'.format(nb_ex_train))


    path_couple_count = f"{path_folder_Result}/table_2d_count"
    np.save(path_couple_count , count_coupleAA)

    return count_AA, nb_AA, count_coupleAA, nb_coupleAA




def freq_for_table_2d(count_AA, nb_AA, count_coupleAA, nb_coupleAA, path_folder_Result):
    """
    Compute and save the frequences of each valid amino acid
    and each valid couple of amino acids.
    """
    # get the list of valid residus
    list_residu = count_AA.keys()

    t = Timer()
    t.start()

    # frequence of each valid amino acid
    freq_AA = {}
    for aa in list_residu:
        if nb_AA != 0:
            freq_AA[aa] = count_AA[aa]/nb_AA
        else: 
            freq_AA[aa] = 0
    
    # frequence of each couple of amino acid
    freq_coupleAA = {}
    for aa_1 in list_residu:
        freq_coupleAA[aa_1] = {}
        for aa_2 in list_residu:
            if nb_coupleAA != 0:
                freq_coupleAA[aa_1][aa_2] = count_coupleAA[aa_1][aa_2]/nb_coupleAA
            else:
                freq_coupleAA[aa_1][aa_2] = 0
    print("")
    t.stop("calcul de table_2d_freq")

    path_freqAA = f"{path_folder_Result}/table_2d_freq"
    np.save(path_freqAA, freq_AA)

    return freq_AA, freq_coupleAA


def table_2d_score(freq_AA, freq_coupleAA, path_folder_Result, scale_factor = 2):
    """
    Compute and save table_2d_score 
    """
    list_residu = freq_AA.keys()

    t = Timer()
    t.start()
    table_2d = {}
    for aa_1 in list_residu:
        table_2d[aa_1] = {}
        for aa_2 in list_residu:
            if freq_coupleAA[aa_1][aa_2] != 0:
                table_2d[aa_1][aa_2] = round(scale_factor * log2(freq_coupleAA[aa_1][aa_2]
                                                            / (freq_AA[aa_1] * freq_AA[aa_2])))
            else:
                table_2d[aa_1][aa_2] = 0
    print("")
    t.stop("calcul de table_2d_score")

    path_matrix = f"{path_folder_Result}/table_2d_score"
    np.save(path_matrix, table_2d) 

    return table_2d


def table_2d_conditional_proba(freq_AA, freq_coupleAA, path_folder_Result):
    """
    Compute and save the matrix of conditional probabilities
    """
    list_residu = freq_AA.keys()

    t = Timer()
    t.start()
    cond_proba = {}
    for aa_1 in list_residu:
        cond_proba[aa_1] = {}
        for aa_2 in list_residu:
            if freq_coupleAA[aa_1][aa_2] != 0:
                cond_proba[aa_1][aa_2] = freq_coupleAA[aa_1][aa_2]/freq_AA[aa_1]
            else:
                cond_proba[aa_1][aa_2] = 0
    print("")
    t.stop("calcul de table_2d_proba")

    path_cond_proba = f"{path_folder_Result}/table_2d_proba"
    np.save(path_cond_proba, cond_proba)

    return cond_proba


def table_2d_heatmap(matrix, path_folder_Result, title, size_annot = 3):
    """
    Save the heatmap of the matrix in path_folder_Result
    """
    #heatmap_matrix = pd.DataFrame(matrix).T.fillna(0)
    heatmap_matrix = np.transpose(pd.DataFrame.from_dict(matrix))
    #cmap = sb.diverging_palette(145, 300, s=60, as_cmap=True)  # test de palette de couleurs
    #cmap = sb.color_palette("vlag", as_cmap=True)
    #heatmap = sb.heatmap(heatmap_matrix, annot = True, annot_kws = {"size": size_annot}, fmt = '.2g', cmap = cmap)
    heatmap = sb.heatmap(heatmap_matrix, annot = True, annot_kws = {"size": size_annot}, fmt = '.2g')    
    plt.yticks(rotation=0) 
    heatmap_figure = heatmap.get_figure()    
    plt.title(title, loc='center', wrap=True)
    #plt.title(title)
    plt.close()
    path_save_fig = f"{path_folder_Result}/{title}.png"
    heatmap_figure.savefig(path_save_fig, dpi=400)



def table_2d_visualisation_transposition(table_2d):
    """
    Visualisation of the matrix
    """
    df_table_2d = np.transpose(pd.DataFrame.from_dict(table_2d))
    print(df_table_2d)

    return df_table_2d


def sum_line_transposition(table_2d):
    """
    To check that the sum of a line is equal to one
    for the conditional probability matrix
    """
    df_table_2d= np.transpose(pd.DataFrame.from_dict(table_2d))
    sum_line = df_table_2d.sum(axis=1)
    print("Somme de chaque ligne :")
    print("")
    print(sum_line)


def table_2d_difference(table_2d, pid_inf_ref):
    """
    Quantify the distance between the table_2d_score computed and a blosum of reference
    """

    list_residu = table_2d.keys()

    # blosum ref importation
    blosum_ref = bl.BLOSUM(pid_inf_ref) 

    # initialisation
    matrix_diff = {}
    difference = 0
    count = 0

    # evaluation of the differences
    for aa1 in list_residu:
        matrix_diff[aa1] = {}
        for aa2 in list_residu:
            matrix_diff[aa1][aa2] = int(table_2d[aa1][aa2] - blosum_ref[aa1 + aa2])
            difference += matrix_diff[aa1][aa2]
            count += 1
    average_difference  = round(difference/count, 2)

    return matrix_diff, pid_inf_ref, average_difference

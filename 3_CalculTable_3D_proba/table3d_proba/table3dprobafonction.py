################################################################################
#                                  Importations                                #    
################################################################################
import numpy as np


import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1]) 

from utils.timer import Timer


################################################################################
#                                  Fonctions                                   #    
################################################################################



# def table_3d_proba_m1(table_3d_count, delay_num, kp_SeqChoice, 
#                       ALPHABET,
#                       path_NeighborRes_m1):
#     """
#     Calcul de la tables 3d proba (delay_num, kp_SeqChoice) de la table 3d count selon la méthode 1
#     Puis sauvegarde de cette table. 
#     """
#     t = Timer()
#     t.start()

#     intra_couple_count = {}
#     for aa_k in ALPHABET:
#         intra_couple_count[aa_k] = {}
#         for aa_c in ALPHABET: 
#             intra_couple_count[aa_k][aa_c] = 0
#             for aa_p in ALPHABET:
#                 intra_couple_count[aa_k][aa_c] += table_3d_count[aa_k][aa_p][aa_c]
#     table_3d_proba_m1 = {}
#     for aa_k in ALPHABET:
#         table_3d_proba_m1[aa_k] = {}
#         for aa_p in ALPHABET:
#             table_3d_proba_m1[aa_k][aa_p] = {}
#             for aa_c in ALPHABET: 
#                 if intra_couple_count[aa_k][aa_c] != 0:
#                     table_3d_proba_m1[aa_k][aa_p][aa_c] = (table_3d_count[aa_k][aa_p][aa_c]) / (intra_couple_count[aa_k][aa_c])  
#                 else:
#                     table_3d_proba_m1[aa_k][aa_p][aa_c] = 0

#     path_table_3d_proba_m1 = f"{path_NeighborRes_m1}/table_3d_proba_M1_({str(delay_num)},{kp_SeqChoice})"
#     np.save(path_table_3d_proba_m1 , table_3d_proba_m1) 
#     path_table_3d_proba_m1  = f"{path_table_3d_proba_m1 }.npy"

#     t.stop("temps de calcul")

#     return table_3d_proba_m1, path_table_3d_proba_m1                 




def table_3d_proba(table_3d_count, delay_num, kp_SeqChoice, 
                    ALPHABET,
                    path_NeighborRes):
    """
    Calcul de la tables 3d proba (delay_num, kp_SeqChoice) de la table 3d count selon la méthode 2
    Puis sauvegarde de cette table. 
    """
    t = Timer()
    t.start()

    intra_couple_count = {}
    for aa_k in ALPHABET:
        intra_couple_count[aa_k] = {}
        for aa_p in ALPHABET: 
            intra_couple_count[aa_k][aa_p] = 0
            for aa_c in ALPHABET:
                intra_couple_count[aa_k][aa_p] += table_3d_count[aa_k][aa_p][aa_c]

    table_3d_proba_m2 = {}
    for aa_k in ALPHABET:
        table_3d_proba_m2[aa_k] = {}
        for aa_p in ALPHABET:
            table_3d_proba_m2[aa_k][aa_p] = {}
            for aa_c in ALPHABET: 
                if intra_couple_count[aa_k][aa_p] != 0:
                    table_3d_proba_m2[aa_k][aa_p][aa_c] = (table_3d_count[aa_k][aa_p][aa_c]) / (intra_couple_count[aa_k][aa_p])  
                else:
                    table_3d_proba_m2[aa_k][aa_p][aa_c] = 0

    path_table_3d_proba_m2 = f"{path_NeighborRes}/table_3d_proba_({str(delay_num)},{kp_SeqChoice})"
    np.save(path_table_3d_proba_m2 , table_3d_proba_m2) 
    path_table_3d_proba_m2  = f"{path_table_3d_proba_m2}.npy"

    t.stop("temps de calcul")

    return table_3d_proba_m2, path_table_3d_proba_m2   




# def sum_line_m1(table_3d_proba_m1, ALPHABET, aa_k, aa_c):
#     """
#     The sum of the conditional probabilities on each line (for aa_k and aa_c fixed) 
#     must be equal to 1 to respect the total probability formula.
#     """
#     sum_line = 0
#     for aa_p in ALPHABET:
#         sum_line += table_3d_proba_m1[aa_k][aa_p][aa_c]
#     print(sum_line)
#     return sum_line


def sum_line(table_3d_proba_m2, ALPHABET, aa_k, aa_p):
    """
    The sum of the conditional probabilities on each line (for aa_k and aa_c fixed) 
    must be equal to 1 to respect the total probability formula.
    """
    sum_line = 0
    for aa_c in ALPHABET:
        sum_line += table_3d_proba_m2[aa_k][aa_p][aa_c]
    print(sum_line)
    return sum_line



def sum_plate(table_3d_proba):
    """
    cond_proba: cube of the conditional probabilities of each valid triplet.

    The sum of the conditional probabilities on each horizontal level of the cube
    must be equal to the length of an edge of the cube to respect the total probability formula.
    """
    for aa_1 in table_3d_proba:
        sum_plateau = 0
        for aa_2 in table_3d_proba[aa_1]:
            for aa_3 in table_3d_proba[aa_1][aa_2]:
                sum_plateau += table_3d_proba[aa_1][aa_2][aa_3]
        print(f"{aa_1}, {sum_plateau}")
"""
TEST: Brier Score
"""

# IMPORTS
import numpy as np
import time

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])



def table_3d_proba_loader(max_relative_distance,
                          origine_destination, left_right,
                          path_table_3d_proba_folder):
    """Load the table_3d_proba of a context

    Args:
        max_relative_distance (float): max position distance of the contextual amino-acid
                                       with the reference
        origine_destination (str): pick the reference as origine 'o' or destination 'd'
        left_right (str): pick left or right as the direction of interest
        path_table_3d_proba_folder (str): path to the folder of all the table_3d_proba

    Returns:
        list: table_3d_proba of interest
    """
    list_table_3d_proba_quarter_window = []

    for i in range(1, max_relative_distance + 1):
        if left_right == "left":
            file_name = f"-{i}_{origine_destination}.npy"
        if left_right == "right":
            file_name = f"{i}_{origine_destination}.npy"
        
        path_cube = f"{path_table_3d_proba_folder}/{file_name}"
        list_table_3d_proba_quarter_window.append(np.load(path_cube, allow_pickle='TRUE').item())

    return list_table_3d_proba_quarter_window
    


def vecteur_from_table_3d_proba(aa_origine, context_str, list_table_3d_proba, alphabet):
    """Retreive vectors from table_3d_proba for a defined context_str
       In case the contextual amino-acid is not in alphabet,
       the corresponding vector is the unit vector

    Args:
        aa_origine (str): origine amino_acid
        context_str (str): contextual amino_acids
        list_table_3d_proba (list): list of the table_3d_proba
        alphabet (list): list of characters of interest

    Returns:
        list: list of vectors
    """
    # print(f"AA ORIGINE: {aa_origine}")
    # print(f"CONTEXT: {context_str}")
    
    len_alphabet = len(alphabet)
    list_vect = [np.array(len_alphabet*[1])]
    for index, aa_voisin in enumerate(context_str):
        if aa_voisin in alphabet:
            vect = []
            for aa in alphabet:
                vect.append(list_table_3d_proba[index][aa_origine][aa][aa_voisin])
            list_vect.append(np.array(vect))
        else:
            list_vect.append(np.array(len_alphabet*[1]))
    # print(f"VECTORS: {list_vect}")
    return list_vect





def unit_brier_naive_bayes(vect, aa_destination, alphabet):
    """Compute the brier score for a vector of prediction
       and the correct prediction aa_destination

    Args:
        vect (list): probability of mutation in each element of alphabet
        aa_destination (str): correct prediction
        alphabet (list): list of characters of interest

    Returns:
        float: brier score
    """
    unit_brier = 0
    for index, aa in enumerate(alphabet):
        unit_brier += (vect[index] - int(aa_destination == aa))**2
    return unit_brier



 

def brier_score(list_example,
                context_ol, context_or,
                context_dl, context_dr,
                alphabet,
                path_table_3d_proba_folder, path_table_2d_proba,
                method):
    """Compute Brier score of one test

    Args:
        list_example (list): _the example tests
        context_ol (int): index position of the contextual amino-acid
                          with respect to origine/left
        context_or (int): index position of the contextual amino-acid
                          with respect to origine/right
        context_dl (int): index position of the contextual amino-acid
                          with respect to destination/left
        context_dr (int): index position of the contextual amino-acid
                          with respect to destination/right
        alphabet (list): list of characters of interest
        path_table_3d_proba_folder (str): path to the table_3d_proba
        path_table_2d_proba (str): path to the table_2d_proba
        method (str): 'multi' considering all the contextual amino-acids of the context
                       OR
                       'uni' considering the farthest contextual amino-acid of the context

    Returns:
        float: Brier Score (mean of elementary tests)
        list: list of the Brier Score of each elementary test
    """
    # context validity check (NOT SURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
    count_0 = 0
    list_context = [context_ol, context_or, context_dl, context_dr]
    for elem in list_context:
        if elem == 0:
            count_0 += 1
    if count_0 < 3:
        print("INVALID CONTEXT, ONE QUARTER CONTEXT MAX SHOULD BE NOT NULL")
        return None

    # INITIALISATIONS
    score_brier_naive_bayes  = 0
    list_unit_score_brier = []
    # n_conservation = 0
    # n_conservation_true = 0
    # n_conservation_false = 0
    # n_substitution = 0
    # n_substitution_true = 0
    # n_substitution_false = 0

    # LOADING: 2D/3D-PROBA
    table_2d_proba = np.load(path_table_2d_proba, allow_pickle='TRUE').item()
    list_cube_quarter_window_ol = table_3d_proba_loader(context_ol, "o", "left",
                                                        path_table_3d_proba_folder)
    list_cube_quarter_window_or = table_3d_proba_loader(context_or, "o", "right",
                                                        path_table_3d_proba_folder)
    list_cube_quarter_window_dl = table_3d_proba_loader(context_dl, "d", "left",
                                                        path_table_3d_proba_folder)
    list_cube_quarter_window_dr = table_3d_proba_loader(context_dr, "d", "right",
                                                        path_table_3d_proba_folder)

    start = time.time()
    for example in list_example:
        total_list_vect = []
        aa_1 = example[0]     # amino-acid ref origine
        aa_2 = example[1]     # amino-acid ref destination
        aa_c_ol = example[2]  # contextual amino-acids: origine_left
        aa_c_or = example[3]  # contextual amino-acids: origine_right
        aa_c_dl = example[4]  # contextual amino-acids: destination_left
        aa_c_dr = example[5]  # contextual amino-acids: destination_right

        # Non-contextual vector (a priori distribution)
        vect_distribution = []
        for aa in alphabet:
            vect_distribution.append(table_2d_proba[aa_1][aa])

        # Contextual vectors
        # print("\nORIGINE-LEFT")
        list_vect_ol = vecteur_from_table_3d_proba(aa_1, aa_c_ol, list_cube_quarter_window_ol, alphabet)
        # print("\nORIGINE-RIGHT")
        list_vect_or = vecteur_from_table_3d_proba(aa_1, aa_c_or, list_cube_quarter_window_or, alphabet)
        # print("\nDESTINATION-LEFT")
        list_vect_dl = vecteur_from_table_3d_proba(aa_1, aa_c_dl, list_cube_quarter_window_dl, alphabet)
        # print("\nDESTINATION-RIGHT")
        list_vect_dr = vecteur_from_table_3d_proba(aa_1, aa_c_dr, list_cube_quarter_window_dr, alphabet)

        # VECTORS CONCATENATION
        if method =="multi":
            total_list_vect = [vect_distribution] + list_vect_ol + list_vect_or + list_vect_dl + list_vect_dr
            # print(f"METHOD: {method}")
            # print(f"ALL THE VECTORS: {total_list_vect}")
        if method =="uni":
            total_list_vect = [vect_distribution, list_vect_ol[-1], list_vect_or[-1], list_vect_dl[-1], list_vect_dr[-1]]
            # print(f"METHOD: {method}")
            # print(f"ALL THE VECTORS: {total_list_vect}")

        # VECTORS ELEMENT WISE PRODUCT
        final_vector = np.prod(np.vstack(total_list_vect), axis=0)
        # print(f"FINAL VECTOR: {final_vector}")

        # VECTOR NORMALIZATION
        final_vector_normalized = final_vector/np.sum(final_vector)
        # print(f"FINAL VECTOR NORMALIZED: {final_vector_normalized}")
        # print(f"COORDINATE SUM: {np.sum(final_vector_normalized)}")

        # BRIER SCORE / EXAMPLE
        score_brier_one_example = unit_brier_naive_bayes(final_vector_normalized, aa_2, alphabet)
        list_unit_score_brier.append(score_brier_one_example)
        score_brier_naive_bayes += score_brier_one_example

        # # PREDICTION CHECK
        # aa_origine = aa_1
        # aa_destination = aa_2
        # final_vector_normalized = list(final_vector_normalized)
        # max_proba = max(final_vector_normalized)
        # index_max = final_vector_normalized.index(max_proba)
        # aa_predicted = alphabet[index_max]

        # if aa_predicted == aa_origine: # prediction of a conservation
        #     n_conservation += 1
        #     if aa_predicted == aa_destination:
        #         n_conservation_true +=1
        #     else:
        #         n_conservation_false +=1
        # else: # prediction of a substitution
        #     n_substitution += 1
        #     if aa_predicted == aa_destination:
        #         n_substitution_true +=1
        #     else:
        #         n_substitution_false +=1



    nb_example = len(list_example)
    if nb_example != 0:
        score_brier_naive_bayes /= nb_example

    # perc_conservation_true = 100*n_conservation_true/nb_example
    # perc_conservation_false = 100*n_conservation_false/nb_example
    # perc_substitution_true = 100*n_substitution_true/nb_example
    # perc_substitution_false = 100*n_substitution_false/nb_example

    # prediction_info = (perc_conservation_true, perc_conservation_false,
    #                    perc_substitution_true, perc_substitution_false)

    end = time.time()
    diff = end - start
    items_per_second = nb_example/diff
    print(f'BRIER SCORE: {score_brier_naive_bayes} | time {diff:.2f}s | {items_per_second:.2f}it/s')

    return score_brier_naive_bayes, list_unit_score_brier, nb_example  #, prediction_info






def brier_score_check_convergence(list_example,
                context_ol, context_or,
                context_dl, context_dr,
                alphabet,
                path_table_3d_proba_folder, path_table_2d_proba,
                method,
                vector_diff):
    """Compute Brier score of one test

    Args:
        list_example (list): _the example tests
        context_ol (int): index position of the contextual amino-acid
                          with respect to origine/left
        context_or (int): index position of the contextual amino-acid
                          with respect to origine/right
        context_dl (int): index position of the contextual amino-acid
                          with respect to destination/left
        context_dr (int): index position of the contextual amino-acid
                          with respect to destination/right
        alphabet (list): list of characters of interest
        path_table_3d_proba_folder (str): path to the table_3d_proba
        path_table_2d_proba (str): path to the table_2d_proba
        method (str): 'multi' considering all the contextual amino-acids of the context
                       OR
                       'uni' considering the farthest contextual amino-acid of the context

    Returns:
        float: Brier Score (mean of elementary tests)
        list: list of the Brier Score of each elementary test
    """
    # context validity check (NOT SURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
    count_0 = 0
    list_context = [context_ol, context_or, context_dl, context_dr]
    for elem in list_context:
        if elem == 0:
            count_0 += 1
    if count_0 < 3:
        print("INVALID CONTEXT, ONE QUARTER CONTEXT MAX SHOULD BE NOT NULL")
        return None

    # INITIALISATIONS
    score_brier_naive_bayes  = 0
    list_unit_score_brier = []
    # n_conservation = 0
    # n_conservation_true = 0
    # n_conservation_false = 0
    # n_substitution = 0
    # n_substitution_true = 0
    # n_substitution_false = 0

    # LOADING: 2D/3D-PROBA
    table_2d_proba = np.load(path_table_2d_proba, allow_pickle='TRUE').item()
    list_cube_quarter_window_ol = table_3d_proba_loader(context_ol, "o", "left",
                                                        path_table_3d_proba_folder)
    list_cube_quarter_window_or = table_3d_proba_loader(context_or, "o", "right",
                                                        path_table_3d_proba_folder)
    list_cube_quarter_window_dl = table_3d_proba_loader(context_dl, "d", "left",
                                                        path_table_3d_proba_folder)
    list_cube_quarter_window_dr = table_3d_proba_loader(context_dr, "d", "right",
                                                        path_table_3d_proba_folder)

    start = time.time()
    for example in list_example:
        total_list_vect = []
        aa_1 = example[0]     # amino-acid ref origine
        aa_2 = example[1]     # amino-acid ref destination
        aa_c_ol = example[2]  # contextual amino-acids: origine_left
        aa_c_or = example[3]  # contextual amino-acids: origine_right
        aa_c_dl = example[4]  # contextual amino-acids: destination_left
        aa_c_dr = example[5]  # contextual amino-acids: destination_right

        # Non-contextual vector (a priori distribution)
        vect_distribution = []
        for aa in alphabet:
            vect_distribution.append(table_2d_proba[aa_1][aa])

        # Contextual vectors
        # print("\nORIGINE-LEFT")
        list_vect_ol = vecteur_from_table_3d_proba(aa_1, aa_c_ol, list_cube_quarter_window_ol, alphabet)
        # print("\nORIGINE-RIGHT")
        list_vect_or = vecteur_from_table_3d_proba(aa_1, aa_c_or, list_cube_quarter_window_or, alphabet)
        # print("\nDESTINATION-LEFT")
        list_vect_dl = vecteur_from_table_3d_proba(aa_1, aa_c_dl, list_cube_quarter_window_dl, alphabet)
        # print("\nDESTINATION-RIGHT")
        list_vect_dr = vecteur_from_table_3d_proba(aa_1, aa_c_dr, list_cube_quarter_window_dr, alphabet)

        # VECTORS CONCATENATION
        if method =="multi":
            total_list_vect = [vect_distribution] + list_vect_ol + list_vect_or + list_vect_dl + list_vect_dr
            # print(f"METHOD: {method}")
            # print(f"ALL THE VECTORS: {total_list_vect}")
        if method =="uni":
            total_list_vect = [vect_distribution, list_vect_ol[-1], list_vect_or[-1], list_vect_dl[-1], list_vect_dr[-1]]
            # print(f"METHOD: {method}")
            # print(f"ALL THE VECTORS: {total_list_vect}")

        # VECTORS ELEMENT WISE PRODUCT
        final_vector = np.prod(np.vstack(total_list_vect), axis=0)
        # print(f"FINAL VECTOR: {final_vector}")

        # VECTOR NORMALIZATION
        final_vector_normalized = final_vector/np.sum(final_vector)
        # print(f"FINAL VECTOR NORMALIZED: {final_vector_normalized}")
        # print(f"COORDINATE SUM: {np.sum(final_vector_normalized)}")

        

        # CHECK CONVERGENCE
        for index in range(len(final_vector)):
            vector_diff += abs(final_vector[index] - vect_distribution[index])



        # BRIER SCORE / EXAMPLE
        score_brier_one_example = unit_brier_naive_bayes(final_vector_normalized, aa_2, alphabet)
        list_unit_score_brier.append(score_brier_one_example)
        score_brier_naive_bayes += score_brier_one_example

        # # PREDICTION CHECK
        # aa_origine = aa_1
        # aa_destination = aa_2
        # final_vector_normalized = list(final_vector_normalized)
        # max_proba = max(final_vector_normalized)
        # index_max = final_vector_normalized.index(max_proba)
        # aa_predicted = alphabet[index_max]

        # if aa_predicted == aa_origine: # prediction of a conservation
        #     n_conservation += 1
        #     if aa_predicted == aa_destination:
        #         n_conservation_true +=1
        #     else:
        #         n_conservation_false +=1
        # else: # prediction of a substitution
        #     n_substitution += 1
        #     if aa_predicted == aa_destination:
        #         n_substitution_true +=1
        #     else:
        #         n_substitution_false +=1



    nb_example = len(list_example)
    if nb_example != 0:
        score_brier_naive_bayes /= nb_example

    # perc_conservation_true = 100*n_conservation_true/nb_example
    # perc_conservation_false = 100*n_conservation_false/nb_example
    # perc_substitution_true = 100*n_substitution_true/nb_example
    # perc_substitution_false = 100*n_substitution_false/nb_example

    # prediction_info = (perc_conservation_true, perc_conservation_false,
    #                    perc_substitution_true, perc_substitution_false)

    end = time.time()
    diff = end - start
    items_per_second = nb_example/diff
    print(f'BRIER SCORE: {score_brier_naive_bayes} | time {diff:.2f}s | {items_per_second:.2f}it/s')

    return score_brier_naive_bayes, list_unit_score_brier, nb_example, vector_diff  #, prediction_info






def brier_score_no_context(list_example,
                           alphabet,
                           path_table_2d_proba):


    # INITIALISATIONS
    score_brier_naive_bayes  = 0
    list_unit_score_brier = []

    # LOADING: 2D_PROBA
    table_2d_proba = np.load(path_table_2d_proba, allow_pickle='TRUE').item()

    start = time.time()
    for example in list_example:
        aa_1 = example[0]     # amino-acid ref origine
        aa_2 = example[1]     # amino-acid ref destination

        # Non-contextual vector (a priori distribution)
        vect_distribution = []
        for aa in alphabet:
            vect_distribution.append(table_2d_proba[aa_1][aa])

        # BRIER SCORE / EXAMPLE
        score_brier_one_example = unit_brier_naive_bayes(vect_distribution, aa_2, alphabet)
        list_unit_score_brier.append(score_brier_one_example)
        score_brier_naive_bayes += score_brier_one_example

    nb_example = len(list_example)
    if nb_example != 0:
        score_brier_naive_bayes /= nb_example

    end = time.time()
    diff = end - start
    items_per_second = nb_example/diff
    print(f'BRIER SCORE: {score_brier_naive_bayes} | time {diff:.2f}s | {items_per_second:.2f}it/s')

    return score_brier_naive_bayes, list_unit_score_brier, nb_example

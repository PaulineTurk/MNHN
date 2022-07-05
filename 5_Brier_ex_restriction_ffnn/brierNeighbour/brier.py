################################################################################
#                                  Importations                                #    
################################################################################
import numpy as np
from tqdm import tqdm
import time

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1]) 
from utils.timer import Timer
import pickle



def table_3d_proba_loader(max_relative_distance, k_or_p, l_or_r, 
                          path_table_3d_proba_folder):
    """
    max_relative_distance: indice le plus lointain dans quart de fenetre
    k_or_p: "origine" (known, k) ou "destination" (predicted, p)
    l_or_r: voisinage à gauche ou à droite
    """
    # initialisation de la liste des cubes pour un quart de fenetre contextuelle
    list_table_3d_proba_quarter_window = []

    for i in range(1, max_relative_distance + 1): 
        if l_or_r == "l":
            path_cube = f"{path_table_3d_proba_folder}/table_3d_proba_(-{i},{k_or_p}).npy" 
            list_table_3d_proba_quarter_window.append(np.load(path_cube, allow_pickle='TRUE').item())

        if l_or_r == "r":  # je sais qu'un else marche mais if pour le moment
            path_cube = f"{path_table_3d_proba_folder}/table_3d_proba_({i},{k_or_p}).npy" 
            list_table_3d_proba_quarter_window.append(np.load(path_cube, allow_pickle='TRUE').item())
        
    #print("list_ctable_3d_proba_quarter_window\n", list_table_3d_proba_quarter_window)

    return list_table_3d_proba_quarter_window
    

def vecteur_from_table_3d_proba(aa_origine, contexte, list_table_3d_proba, ALPHABET):
    """
    aa_origine: acide amnié d'origine'
    contexte: liste des caractères d'une fraction de voisinage local (lecture en s'éloigant de l'acide aminé de référence)
    list_cube: list des cubes pré-calculés et pré-chargés nécessaires pour ce contexte
    """
    len_ALPHABET = len(ALPHABET)
    list_vect = [np.array(len_ALPHABET*[1])]
    for index, aa_voisin in enumerate(contexte):
        if aa_voisin in ALPHABET:
            vect = []
            for aa in ALPHABET:
                vect.append(list_table_3d_proba[index][aa_origine][aa][aa_voisin])
            list_vect.append(np.array(vect))
        else:
            list_vect.append(np.array(len_ALPHABET*[1]))

    #print("list_vect: \n",list_vect)

    return list_vect





def unit_brier_naive_bayes(vect, aa_destination, ALPHABET):
    """
    vect: vecteur de distribution de probabilité pré-calculé pour chaque exemple selon son contexte
    aa_2: acide aminé de destination
    """
    unit_brier = 0
    for index, aa in enumerate(ALPHABET):
        unit_brier += (vect[index] - int(aa_destination == aa))**2
    #print("unit_brier:", unit_brier)
    return unit_brier


def brier_score_multi(list_example, 
                    context_kl, context_kr, context_pl, context_pr, 
                    ALPHABET,
                    path_table_3d_proba_folder, path_table_2d_proba):
    """
    list_example; liste des examples selectionnés à tester
    context_kl, context_kr, context_pl, context_pr: indice max de chaque quart de fenetre
    path_table_3d_proba_folder: chemin du dossier contenant les table_3d_proba
    path_table_2d_proba: chemin de la table 2d probabiliste

    path_save_file: ou enregistrer la liste des scores de brier par exemple
    """
    # verification contexte valide
    # seuls les contexts par quart de fenetre sont testés (ou le cas non contextuel)
    count_0 = 0
    list_context = [context_kl, context_kr, context_pl, context_pr]
    for elem in list_context:
        if elem == 0:
            count_0 += 1
    if count_0 < 3:
        #print("count_0:\n", count_0)
        #print("contexte invalide")
        print("invalid context, une unique valeur peut etre non nulle")
        return None
   

    # initialisation du score de Brier Bayes Naif
    score_brier_naive_bayes  = 0

    # initialisation de la liste des scores de brier par exemple
    list_unit_score_brier = []

    # chargement des table 2d et 3d de proba
    table_2d_proba = np.load(path_table_2d_proba, allow_pickle='TRUE').item()
    list_cube_quarter_window_kl = table_3d_proba_loader(context_kl, "origine", "l", path_table_3d_proba_folder)
    list_cube_quarter_window_kr = table_3d_proba_loader(context_kr, "origine", "r", path_table_3d_proba_folder)
    list_cube_quarter_window_pl = table_3d_proba_loader(context_pl, "destination", "l", path_table_3d_proba_folder)
    list_cube_quarter_window_pr = table_3d_proba_loader(context_pr, "destination", "r", path_table_3d_proba_folder)

    #print(list_example)

    for example in list_example:
        #print("\nexample:\n", example)
        total_list_vect = []

        aa_1 = example[0]     # origine
        aa_2 = example[1]     # destination
        aa_c_kl = example[2]  # voisin/s fenetre kl
        aa_c_kr = example[3]  # voisin/s fenetre kr
        aa_c_pl = example[4]  # voisin/s fenetre pl
        aa_c_pr = example[5]  # voisin/s fenetre pl

        # initialisation du vecteur de distribution de probabilité de mutation 
        # de aa_1 en chacun des 20 aa standards
        vect_distribution = []
        for aa in ALPHABET:
            vect_distribution.append(table_2d_proba[aa_1][aa])    # vecteur a priori (non-contextuel)
        
        total_list_vect.append(np.array(vect_distribution))
        #print("init total list with blosum proba:\n", total_list_vect)

        list_vect_kl = vecteur_from_table_3d_proba(aa_1, aa_c_kl, list_cube_quarter_window_kl, ALPHABET)
        #print("list_vect_kl:\n", list_vect_kl)
        list_vect_kr = vecteur_from_table_3d_proba(aa_1, aa_c_kr, list_cube_quarter_window_kr, ALPHABET)
        #print("list_vect_kr:\n", list_vect_kr)
        list_vect_pl = vecteur_from_table_3d_proba(aa_1, aa_c_pl, list_cube_quarter_window_pl, ALPHABET)
        #print("list_vect_pl:\n", list_vect_pl)
        list_vect_pr = vecteur_from_table_3d_proba(aa_1, aa_c_pr, list_cube_quarter_window_pr, ALPHABET)        
        #print("list_vect_pr:\n", list_vect_pr)
        
        #print("list_vect 2d:\n", total_list_vect)
        #print("ok")
        total_list_vect = total_list_vect + list_vect_kl + list_vect_kr + list_vect_pl + list_vect_pr   # vect of arrays
        #print("total_list_vect:\n", total_list_vect)

        final_vector = np.prod(np.vstack(total_list_vect), axis=0) 
        final_vector_normalized = final_vector/np.sum(final_vector)

        #print("final_vector_normalized\n", final_vector_normalized)
        #print("sum vector normalized:", np.sum(final_vector_normalized))

        score_brier_one_example = unit_brier_naive_bayes(final_vector_normalized, aa_2, ALPHABET)

        list_unit_score_brier.append(score_brier_one_example)
        score_brier_naive_bayes += score_brier_one_example

    nb_example = len(list_example)
    if nb_example != 0:
        score_brier_naive_bayes /= nb_example
    print("nombre d'exemples tirés:", nb_example)
    print("score_brier_naive_bayes   :", score_brier_naive_bayes)


    return score_brier_naive_bayes, list_unit_score_brier








def brier_score_uni(list_example, 
                    context_kl, context_kr, context_pl, context_pr, 
                    ALPHABET,
                    path_table_3d_proba_folder, path_table_2d_proba):
    """
    list_example; liste des examples selectionnés à tester
    context_kl, context_kr, context_pl, context_pr: indice max de chaque quart de fenetre
    path_table_3d_proba_folder: chemin du dossier contenant les table_3d_proba
    path_table_2d_proba: chemin de la table 2d probabiliste

    path_save_file: ou enregistrer la liste des scores de brier par exemple
    """
    # verification contexte valide
    # seuls les contexts par quart de fenetre sont testés (ou le cas non contextuel)
    count_0 = 0
    list_context = [context_kl, context_kr, context_pl, context_pr]
    for elem in list_context:
        if elem == 0:
            count_0 += 1
    if count_0 < 3:
        #print("count_0:\n", count_0)
        #print("contexte invalide")
        print("invalid context, une unique valeur peut etre non nulle")
        return None
   

    # initialisation du score de Brier Bayes Naif
    score_brier_naive_bayes  = 0

    # initialisation de la liste des scores de brier par exemple
    list_unit_score_brier = []

    # chargement des table 2d et 3d de proba si le voisinage est étudié
    table_2d_proba = np.load(path_table_2d_proba, allow_pickle='TRUE').item()
    list_cube_quarter_window_kl = table_3d_proba_loader(context_kl, "origine", "l", path_table_3d_proba_folder)
    list_cube_quarter_window_kr = table_3d_proba_loader(context_kr, "origine", "r", path_table_3d_proba_folder)
    list_cube_quarter_window_pl = table_3d_proba_loader(context_pl, "destination", "l", path_table_3d_proba_folder)
    list_cube_quarter_window_pr = table_3d_proba_loader(context_pr, "destination", "r", path_table_3d_proba_folder)

    #print(list_example)

    for example in list_example:
        #print("\nexample:\n", example)
        total_list_vect = []

        aa_1 = example[0]     # origine
        aa_2 = example[1]     # destination
        aa_c_kl = example[2]  # voisin/s fenetre kl
        aa_c_kr = example[3]  # voisin/s fenetre kr
        aa_c_pl = example[4]  # voisin/s fenetre pl
        aa_c_pr = example[5]  # voisin/s fenetre pl

        # initialisation du vecteur de distribution de probabilité de mutation 
        # de aa_1 en chacun des 20 aa standards
        vect_distribution = []
        for aa in ALPHABET:
            vect_distribution.append(table_2d_proba[aa_1][aa])    # vecteur a priori (non-contextuel)
        
        total_list_vect.append(np.array(vect_distribution))
        #print("init total list with blosum proba:\n", total_list_vect)

        list_vect_kl = vecteur_from_table_3d_proba(aa_1, aa_c_kl, list_cube_quarter_window_kl, ALPHABET)
        #print("list_vect_kl:\n", list_vect_kl)
        list_vect_kr = vecteur_from_table_3d_proba(aa_1, aa_c_kr, list_cube_quarter_window_kr, ALPHABET)
        #print("list_vect_kr:\n", list_vect_kr)
        list_vect_pl = vecteur_from_table_3d_proba(aa_1, aa_c_pl, list_cube_quarter_window_pl, ALPHABET)
        #print("list_vect_pl:\n", list_vect_pl)
        list_vect_pr = vecteur_from_table_3d_proba(aa_1, aa_c_pr, list_cube_quarter_window_pr, ALPHABET)        
        #print("list_vect_pr:\n", list_vect_pr)
        
        #print("list_vect 2d:\n", total_list_vect)


        # seule partie qui change entre uni et multi
        total_list_vect = [total_list_vect[-1], list_vect_kl[-1], list_vect_kr[-1], list_vect_pl[-1], list_vect_pr[-1]]   # list of arrays
        #print("total_list_vect:\n", total_list_vect)

        final_vector = np.prod(np.vstack(total_list_vect), axis=0) # pas encore normalisé
        final_vector_normalized = final_vector/np.sum(final_vector)

        #print("final_vector_normalized\n", final_vector_normalized)
        #print("sum vector normalized:", np.sum(final_vector_normalized))

        score_brier_one_example = unit_brier_naive_bayes(final_vector_normalized, aa_2, ALPHABET)

        list_unit_score_brier.append(score_brier_one_example)
        score_brier_naive_bayes += score_brier_one_example

    nb_example = len(list_example)
    if nb_example != 0:
        score_brier_naive_bayes /= nb_example
    print("nombre d'exemples tirés:", nb_example)
    print("score_brier_naive_bayes   :", score_brier_naive_bayes)


    return score_brier_naive_bayes, list_unit_score_brier







def estimationError(list_position, list_brier_score, nb_position, nb_test):

    list_position_err       = []
    list_brier_score_mean   = []
    list_brier_score_std    = []

    i = 0
    while i in range(nb_position):
        list_position_err.append(list_position[i])
        estimation_brier = np.array(list_brier_score[i: i+nb_test])
        list_brier_score_mean.append(np.mean(estimation_brier))
        list_brier_score_std.append(np.std(estimation_brier))
        i += nb_test
    
    return np.array(list_position_err), np.array(list_brier_score_mean), np.array(list_brier_score_std)


def writeList(output_list, path_save_file):
    with open(path_save_file, 'wb') as temp:
        pickle.dump(output_list, temp)


def readList(path_save_file):
    with open (path_save_file, 'rb') as temp:
        items = pickle.load(temp)
    return items




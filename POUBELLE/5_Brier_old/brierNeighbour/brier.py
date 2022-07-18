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


def table_3d_proba_loader_1(max_relative_distance, k_or_p, l_or_r, path_table_3d_proba_folder, method):
    """
    max_relative_distance: indice le plus lointain dans quart de fenetre (entier positif)
    k_or_p: séquence d'origine (k) ou de destination (p)
    l_or_r: voisinage à gauche ou à droite
    method = M1 ou M2

    _1: Ne selectionne que la table 3d de l'acide aminé contextuel défini par les paramètres de la fonction
    """
    if l_or_r == "l":
        path_table_3d = f"{path_table_3d_proba_folder}/table_3d_proba_{method}_(-{max_relative_distance},{k_or_p}).npy" 
        table_3d = np.load(path_table_3d, allow_pickle='TRUE').item()

    if l_or_r == "r":  # je sais qu'un else marche mais if pour le moment
        path_table_3d = f"{path_table_3d_proba_folder}/table_3d_proba_{method}_({max_relative_distance},{k_or_p}).npy" 
        table_3d = np.load(path_table_3d, allow_pickle='TRUE').item()

    return table_3d
    

def vecteur_from_table_3d_proba_M1_1(aa_1, context, table_3d, 
                                list_inclusion, 
                                nb_ex_valide, nb_ex_invalide,
                                nb_ex_bord, nb_ex_char_exclus):
    """
    aa_1: acide amnié de départ
    context: acides aminés voisins avec context[-1] l'acide aminé contextuel
    table_3d: table 3d associé à la position de l'acide aminé contextuel observé
    nb_ex_valide: l'ex est valide si l'acide aminé contextuel est dans la list d'inclusion
    nb_ex_invalide : sinon il est invalide

    # il y a 2 catégories d'invalidité de l'exemple:
    nb_ex_bord: l'aa contextuel n'existe pas car en position hors de l'alignement symbole *
    nb_ex_char_exclus: le charactère en position de l'aa_c d'interet n'est pas dans la liste
                       d'inclusion et n'est pas le symbole *
    """

    aa_c = context[-1]
    vect = []

    if aa_c in list_inclusion:
        nb_ex_valide += 1
        for aa in list_inclusion:
            vect.append(table_3d[aa_1][aa][aa_c])
    
    else:
        nb_ex_invalide += 1
        if aa_c == "*": # position hors longueur de l'alignement
            nb_ex_bord += 1
        else:
            nb_ex_char_exclus += 1 # position existe 
                                   # mais le caractère n'est pas dans la liste d'inclusion
    #print(vect, nb_ex_valide, nb_ex_invalide, nb_ex_bord, nb_ex_char_exclus)
    return vect, nb_ex_valide, nb_ex_invalide, nb_ex_bord, nb_ex_char_exclus


def vecteur_from_table_3d_proba_M2_1(aa_1, context, table_3d, 
                                list_inclusion, 
                                nb_ex_valide, nb_ex_invalide,
                                nb_ex_bord, nb_ex_char_exclus):
    aa_c = context[-1]
    vect = []

    if aa_c in list_inclusion:
        nb_ex_valide += 1
        for aa in list_inclusion:
            vect.append(table_3d[aa_1][aa][aa_c])
    
    else:
        nb_ex_invalide += 1
        if aa_c == "*": # position hors longueur de l'alignement
            nb_ex_bord += 1
        else:
            nb_ex_char_exclus += 1 # position existe 
                                   # mais le caractère n'est pas dans la liste d'inclusion
    #print(vect, nb_ex_valide, nb_ex_invalide, nb_ex_bord, nb_ex_char_exclus)
    return vect, nb_ex_valide, nb_ex_invalide, nb_ex_bord, nb_ex_char_exclus


def unit_brier_naive_bayes(vect, aa_2, list_inclusion):
    """
    vect: vecteur de distribution de probabilité pré-calculé pour chaque exemple selon son contexte
    aa_2: acide aminé de destination
    """
    unit_brier = 0
    for index, aa in enumerate(list_inclusion):
        unit_brier += (vect[index] - int(aa_2 == aa))**2
    #print("unit_brier:", unit_brier)
    return unit_brier


def brier_score_m1_1(list_example, 
                     context_kl, context_kr, context_pl, context_pr, 
                     list_inclusion,
                     path_table_3d_proba_folder, path_table_2d_proba,
                     method):
    """
    list_example; liste des examples selectionnés à tester
    context_kl, context_kr, context_pl, context_pr: indice max de chaque quart de fenetre
    path_table_3d_proba_folder: chemin du dossier contenant les table_3d_proba
    path_table_2d_proba: chemin de la table 2d probabiliste

    path_save_file: ou enregistrer la liste des scores de brier par exemple
    """
    # verification contexte valide:
    count_0 = 0
    list_context = [context_kl, context_kr, context_pl, context_pr]
    for index, position in enumerate(list_context):
        if position == 0:
            count_0 += 1
        else: 
            position_voisin = index # position_voisin: 0, 1, 2 ou 3 selon seq o/d et coté g/d
            #print("position_voisin:",position_voisin)
            decalage_position = position  # 1, 2, 3, 4 , 5, 6, ... selon la profondeur de tables 3d calculés
            #print("decalage_position:", decalage_position)
    if count_0 < 3:
        #print("count_0:\n", count_0)
        #print("contexte invalide")
        print("invalid context, une unique valeur peut etre non nulle")
        return None
   

    # initialisation du score de Brier
    score_brier       = 0

    # initialisation de la liste des scores de brier par exemple
    list_unit_score_brier = []

    # initialisation du compte des exemples de test valides/invalides
    nb_ex_valide      = 0
    nb_ex_invalide    = 0
    nb_ex_bord        = 0
    nb_ex_char_exclus = 0
    nb_ex_evaluated   = 0   # vérifier qu'il est égal à nb_ex_valide


    # Chargement des tables 2d/3d    
    if count_0 == 4:  # cas non contextuel, chargement de la table 2d
        table_2d = np.load(path_table_2d_proba, allow_pickle='TRUE').item()
        #print("table_2d:\n", table_2d)

    else: # cas contextuel, chargement de la table 3d d'intéret
        if position_voisin == 0:
            kp_SeqChoice, sens = "k", "l"
        if position_voisin == 1:
            kp_SeqChoice, sens = "k", "r"
        if position_voisin == 2:
            kp_SeqChoice, sens = "p", "l"
        if position_voisin == 3:
            kp_SeqChoice, sens = "p", "r"

        table_3d = table_3d_proba_loader_1(decalage_position, kp_SeqChoice, sens, path_table_3d_proba_folder, method)
        #print("table_3d:\n", table_3d)
    #print(list_example)

    print("")
    start = time.time()
    len_list_example = len(list_example)
    for index_exemple in tqdm(range(len_list_example), desc = "Calcul du score de Brier M1.1",
                                    ncols= 100, mininterval=60):
        example = list_example[index_exemple]
        #print("example:", example)
        aa_1 = example[0]
        aa_2 = example[1]
        aa_c_kl = example[2]
        aa_c_kr = example[3]
        aa_c_pl = example[4]
        aa_c_pr = example[5]
        #print("exemple test:\n", example)

        if count_0 == 4:
            nb_ex_valide += 1 # tous les ex non contextels sont valides par construction
            vect = [] # initialisation du vecteur de proba associé à l'exemple test
            for aa in list_inclusion:
                vect.append(table_2d[aa_1][aa])
            #print("vecteur 2d:\n", vect)

        else:
            if position_voisin == 0:
                context = aa_c_kl
            if position_voisin == 1:
                context = aa_c_kr
            if position_voisin == 2:
                context = aa_c_pl
            if position_voisin == 3:
                context = aa_c_pr

            #print("aa_c:", aa_c)
            vect, nb_ex_valide, nb_ex_invalide, nb_ex_bord, nb_ex_char_exclus = vecteur_from_table_3d_proba_M1_1(aa_1, context, table_3d, 
                                                                                                            list_inclusion,
                                                                                                            nb_ex_valide, nb_ex_invalide,
                                                                                                            nb_ex_bord, nb_ex_char_exclus)   
        #print("vect:", vect)
        if vect != []:
            nb_ex_evaluated += 1 # doit etre égal au nombre d'ex valides
            unit_score_brier = unit_brier_naive_bayes(vect, aa_2, list_inclusion) 
            score_brier += unit_score_brier 
            list_unit_score_brier.append(unit_score_brier)

    # calcul du score de Brier
    if nb_ex_evaluated != 0:
        score_brier /= nb_ex_evaluated
        print("score_brier :", score_brier)
        print("nb d'exemples utilisés pour le calcul du score de Brier:", nb_ex_evaluated)
    else:
        print("nb d'exemples utilisés pour le calcul du score de Brier:", nb_ex_evaluated)
        print("score de Brier non évaluable")


    # bilan du comptage des exemples
    print("\nBilan de la valdité des exemples")
    nb_example = len(list_example)
    if nb_example != 0:
        print(f"nb ex valid                (% total ex): {'{:_}'.format(nb_ex_valide)} ({round(100*nb_ex_valide/nb_example, 2)} %)")   
        print(f"nb ex invalid              (% total ex): {'{:_}'.format(nb_ex_invalide)} ({round(100*nb_ex_invalide/nb_example, 2)} %)")  

    if nb_ex_invalide != 0:
        print(f"nb ex nb_ex_char_exclus (% ex invalide): {'{:_}'.format(nb_ex_char_exclus)} ({round(100*nb_ex_char_exclus/nb_ex_invalide, 2)} %)") 
        print(f"nb ex nb_ex_bord        (% ex invalide): {'{:_}'.format(nb_ex_bord)} ({round(100*nb_ex_bord/nb_ex_invalide, 2)} %)") 

    end = time.time()
    diff = end-start
    items_per_second = nb_ex_evaluated/diff
    print(f'Calcul du score de Brier: {diff:.2f} s | {items_per_second:.2f} it/s')

    return score_brier, list_unit_score_brier



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






def brier_score_m2_1(list_example, 
                        context_kl, context_kr, context_pl, context_pr, 
                        list_inclusion,
                        path_table_3d_proba_folder, path_table_2d_proba,
                        method):
    """
    list_example; liste des examples selectionnés à tester
    context_kl, context_kr, context_pl, context_pr: indice max de chaque quart de fenetre
    path_table_3d_proba_folder: chemin du dossier contenant les table_3d_proba
    path_table_2d_proba: chemin de la table 2d probabiliste

    path_save_file: ou enregistrer la liste des scores de brier par exemple
    """
    # verification contexte valide:
    count_0 = 0
    list_context = [context_kl, context_kr, context_pl, context_pr]
    for index, position in enumerate(list_context):
        if position == 0:
            count_0 += 1
        else: 
            position_voisin = index # position_voisin: 0, 1, 2 ou 3 selon seq o/d et coté g/d
            #print("position_voisin:",position_voisin)
            decalage_position = position  # 1, 2, 3, 4 , 5, 6, ... selon la profondeur de tables 3d calculés
            #print("decalage_position:", decalage_position)
    if count_0 < 3:
        #print("count_0:\n", count_0)
        #print("contexte invalide")
        print("invalid context, une unique valeur peut etre non nulle")
        return None
   

    # initialisation du score de Brier
    score_brier       = 0

    # initialisation de la liste des scores de brier par exemple
    list_unit_score_brier = []

    # initialisation du compte des exemples de test valides/invalides
    nb_ex_valide      = 0
    nb_ex_invalide    = 0
    nb_ex_bord        = 0
    nb_ex_char_exclus = 0
    nb_ex_evaluated   = 0   # vérifier qu'il est égal à nb_ex_valide


    # Chargement de la table 2d de proba 
    table_2d = np.load(path_table_2d_proba, allow_pickle='TRUE').item()
        #print("table_2d:\n", table_2d)

    # chargement de la table 3d de proba si le contexte est différent de (0,0,0,0)
    if  list_context != [0, 0, 0, 0]:
        if position_voisin == 0:
            kp_SeqChoice, sens = "k", "l"
        if position_voisin == 1:
            kp_SeqChoice, sens = "k", "r"
        if position_voisin == 2:
            kp_SeqChoice, sens = "p", "l"
        if position_voisin == 3:
            kp_SeqChoice, sens = "p", "r"

        table_3d = table_3d_proba_loader_1(decalage_position, kp_SeqChoice, sens, path_table_3d_proba_folder, method)
        #print("table_3d:\n", table_3d)
    #print(list_example)

    print("")
    start = time.time()
    len_list_example = len(list_example)
    for index_exemple in tqdm(range(len_list_example), desc = "Calcul du score de Brier M1.1",
                                    ncols= 100, mininterval=60):
        example = list_example[index_exemple]
        #print("example:", example)
        aa_1 = example[0]
        aa_2 = example[1]
        aa_c_kl = example[2]
        aa_c_kr = example[3]
        aa_c_pl = example[4]
        aa_c_pr = example[5]
        #print("exemple test:\n", example)

        if count_0 == 4:
            nb_ex_valide += 1 # tous les ex non contextels sont valides par construction
            vect = [] # initialisation du vecteur de proba associé à l'exemple test
            for aa in list_inclusion:
                vect.append(table_2d[aa_1][aa])
            #print("vecteur 2d:\n", vect)

        else:
            if position_voisin == 0:
                context = aa_c_kl
            if position_voisin == 1:
                context = aa_c_kr
            if position_voisin == 2:
                context = aa_c_pl
            if position_voisin == 3:
                context = aa_c_pr

            #print("aa_c:", aa_c)
            vect, nb_ex_valide, nb_ex_invalide, nb_ex_bord, nb_ex_char_exclus = vecteur_from_table_3d_proba(aa_1, context, table_3d, 
                                                                                                            list_inclusion,
                                                                                                            nb_ex_valide, nb_ex_invalide,
                                                                                                            nb_ex_bord, nb_ex_char_exclus)   
        #print("vect:", vect)
        if vect != []:
            nb_ex_evaluated += 1 # doit etre égal au nombre d'ex valides
            unit_score_brier = unit_brier_naive_bayes(vect, aa_2, list_inclusion) 
            score_brier += unit_score_brier 
            list_unit_score_brier.append(unit_score_brier)

    # calcul du score de Brier
    if nb_ex_evaluated != 0:
        score_brier /= nb_ex_evaluated
        print("score_brier :", score_brier)
        print("nb d'exemples utilisés pour le calcul du score de Brier:", nb_ex_evaluated)
    else:
        print("nb d'exemples utilisés pour le calcul du score de Brier:", nb_ex_evaluated)
        print("score de Brier non évaluable")


    # bilan du comptage des exemples
    print("\nBilan de la valdité des exemples")
    nb_example = len(list_example)
    if nb_example != 0:
        print(f"nb ex valid                (% total ex): {'{:_}'.format(nb_ex_valide)} ({round(100*nb_ex_valide/nb_example, 2)} %)")   
        print(f"nb ex invalid              (% total ex): {'{:_}'.format(nb_ex_invalide)} ({round(100*nb_ex_invalide/nb_example, 2)} %)")  

    if nb_ex_invalide != 0:
        print(f"nb ex nb_ex_char_exclus (% ex invalide): {'{:_}'.format(nb_ex_char_exclus)} ({round(100*nb_ex_char_exclus/nb_ex_invalide, 2)} %)") 
        print(f"nb ex nb_ex_bord        (% ex invalide): {'{:_}'.format(nb_ex_bord)} ({round(100*nb_ex_bord/nb_ex_invalide, 2)} %)") 

    end = time.time()
    diff = end-start
    items_per_second = nb_ex_evaluated/diff
    print(f'Calcul du score de Brier: {diff:.2f} s | {items_per_second:.2f} it/s')

    return score_brier, list_unit_score_brier
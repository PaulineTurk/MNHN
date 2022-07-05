################################################################################
#                                  Importations                                #    
################################################################################
from pathlib import Path
import numpy as np
import pandas as pd
import random
import math
import os
from tqdm import tqdm
import time



import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1]) 

import utils.fastaReader as fastaReader 
import utils.folder as folder
from utils.timer import Timer


def one_seed_selection(accession_num, seed, 
                       pid, pid_inf, pid_sup, 
                       nb_non_info_seed, dico_seed, nb_valid_aa_couple_global, 
                       ALPHABET):
    """
    accession_num: pour retrouver le seed dans le dico
    """

    len_align = len(seed[0][1])   # à calculer qu'une seule fois car une longeur d'alignement dans un seed
    nb_seq = len(seed)
    #print(f"seed accession num, nb seq, len alignment: {accession_num}, {len_align}, {nb_seq}")

    dico_seq = {}
    nb_valid_aa_couple_seed = 0

    for i in range(nb_seq):
        name_A, seq_A = seed[i]   
        for j in range(i+1, nb_seq):
            name_B, seq_B = seed[j]

            if pid[name_A][name_B] >= pid_inf and pid[name_A][name_B] < pid_sup:
                # initialisation de la liste des positions valides hors context
                list_valid_position = []
                #for index, aa_A in enumerate(seq_A):  

                for index in range(len_align):   
                    if seq_A[index] in ALPHABET:
                        if seq_B[index] in ALPHABET:
                            list_valid_position.append(index)

                
                nb_valid_position = len(list_valid_position)
                        

                if nb_valid_position != 0:    # ajouter le couple name_A, name_B au dico_seq en cours de remplissage 
                                              # ssi il contient au moins une séquence non contextuelle valide

                    dico_seq[(name_A, name_B)] = {}
                    dico_seq[(name_A, name_B)]["nbValidPosition"] = nb_valid_position
                    dico_seq[(name_A, name_B)]["SeqAB"] = (seq_A, seq_B)
                    dico_seq[(name_A, name_B)]["validPosition"] = list_valid_position   # liste des positions valides
                                          
                nb_valid_aa_couple_seed += nb_valid_position  # compté sur les paires de séquences possibles et non les couples


    if nb_valid_aa_couple_seed == 0: # seed inutil
        nb_non_info_seed += 1    

    else:
        dico_seed[accession_num] = {}
        dico_seed[accession_num]["nbValidPairAA"] = nb_valid_aa_couple_seed
        dico_seed[accession_num]["lenAlign"] = len_align 

    nb_valid_aa_couple_global += nb_valid_aa_couple_seed

    return dico_seed, nb_non_info_seed, nb_valid_aa_couple_global, dico_seq


def multi_seeds_selection(path_folder_seed, path_folder_pid, 
                          pid_inf, pid_sup, ALPHABET, 
                          path_folder_dico_seq, path_file_dico_seed):
    """
    path_folder_dico_seq: dossier créé s'il n'existe pas, effacé puis recréé sinon.
                          Les dictionnaires de seed y sont enregistrés.
    path_file_dico_seed: dossier créé s'il n'existe pas, effacé puis recréé sinon.
                         Le dictionnaire sur les seeds y est enregistré 
                         ainsi que sa version normalisée par la somme de tous les comptes.
    """

    # initialisation des "compteurs" cumulés sur tous les seeds
    dico_seed = {}   # len(dico_seed) = nb_seed_valid
    nb_non_info_seed = 0  # pour vérifier ce compte + len(dico_seed) = len(path_folder_seed)
    nb_valid_aa_couple_global = 0

    # creation of path_folder_dico_seq
    path_folder_dico_seq = folder.creat_folder(path_folder_dico_seq)

    # creation of path_folder_dico_seed
    path_file_dico_seed = folder.creat_folder(path_file_dico_seed)

    path, dirs, files = next(os.walk(path_folder_seed))
    nb_files = len(files)

    # liste des PosixPath des alignements d'apprentissage
    files = [x for x in Path(path_folder_seed).iterdir()]
    nb_files = len(files)
    print("Nombre d'alignements multiples: ", nb_files)
    start = time.time()
    for file_counter in tqdm(range(nb_files), desc='Selection des alignements valides',
                                   ncols= 100, mininterval=60):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        pid = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()
        seed = fastaReader.read_multi_fasta(file)

        dico_seed, nb_non_info_seed, nb_valid_aa_couple_global, dico_seq = one_seed_selection(accession_num, seed, 
                                                                                              pid, pid_inf, pid_sup, 
                                                                                              nb_non_info_seed, dico_seed, nb_valid_aa_couple_global, 
                                                                                              ALPHABET)
        
        #print(f"\n{accession_num}")
        #print(np.transpose(pd.DataFrame.from_dict(dico_seq)))

        # sauvegarde des dico_seq
        path_file_dico_seq = f"{path_folder_dico_seq}/{accession_num}.seq"
        np.save(path_file_dico_seq, dico_seq) 

    # sauvegarde du dico_seed
    path_save_file_dico_seed = f"{path_file_dico_seed}/seed"
    np.save(path_save_file_dico_seed, dico_seed) 

    # normalisation de dico_seed
    dico_seed_normalised = {}
    for key in dico_seed:
        dico_seed_normalised[key] = {}
        dico_seed_normalised[key]["proportionValidPosition"] = dico_seed[key]["nbValidPairAA"]/nb_valid_aa_couple_global
        dico_seed_normalised[key]["lenAlign"] = dico_seed[key]["lenAlign"] 
    
    #print("dico_seed_normalised")
    #print(np.transpose(pd.DataFrame.from_dict(dico_seed_normalised)))
        
    # sauvegarde du dico_seed_normalised
    path_save_file_dico_seed_normalised = f"{path_file_dico_seed}/seed_normalised"
    np.save(path_save_file_dico_seed_normalised, dico_seed_normalised) 

    print(f"\nnb de seed valides, invalides: {len(dico_seed)}, {nb_non_info_seed}")
    end = time.time()
    diff = end - start
    items_per_second = nb_files/diff
    print(f'Selection des alignements valides: {diff:.2f}s | {items_per_second:.2f}it/s')

    

# def example_number_per_seed(path_dico_seed_normalised, nb_exemple_test, path_dico_exemple):
#     """
#     déterminer pseudo-aléatoirement le nombre d'exemple valide à piocher de chaque seed. 

#     path_dico_seed: pré-calculé (à loader une fois)
#     nb_tes_exemple: ordre de grandeur du nombre de couples valides à piocher au total de Pfam test.

#     path_dico_exemple: chemin au fichier ou le file des exemples à prendre de chaque seed est enregistré
#     """
#     t = Timer()
#     t.start()    

#     dico_seed_normalised = np.load(path_dico_seed_normalised, allow_pickle='TRUE').item()

#     # conversion des proportions associées à chaque seed en nombre d'exemples à prendre de chaque seed

#     ## multiplication de chaque proportion par le nombre d'exemples souhaité
#     dico_exemple = {}
#     for key in dico_seed_normalised:
#         dico_exemple[key] = {}
#         nb_ex_estimation = dico_seed_normalised[key]["proportionValidPosition"] * nb_exemple_test

#         # partie entière inf + int(random.uniform(0, 1) < partie décimale)      
#         decimal_part, integer_part = math.modf(nb_ex_estimation)
#         proba = random.uniform(0, 1)
#         nb_ex_exact = int(integer_part + int(proba < decimal_part))   # généralisation en prévision de l'envoie d'exemple par batch 
#                                                                      # au réseau de neurones (respecter l'ordre de grandeur du batch demandé ...) 
#         dico_exemple[key]["nbExample"] = nb_ex_exact
#         dico_exemple[key]["lenAlign"] = dico_seed_normalised[key]["lenAlign"]

#     path_dico_exemple = f"{path_dico_exemple}/exemple_seed"
#     np.save(path_dico_exemple, dico_exemple) 

#     #print(np.transpose(pd.DataFrame.from_dict(dico_exemple)))
#     print("")
#     t.stop("Selection du nombre d'exemples par seed")




# def random_example_selection(list_example, dico_seq, 
#                              context_kl, context_kr, context_pl, context_pr, 
#                              nb_ex_test):   # rq. ex_bord inclus dans ex_test
#     """
#     A ce stade on suppose qu on a selectionné un seed pour y piocher ce qu'on veut et on met cette fonction
#     dans une boucle autant de fois qu'on doit prendre d'exemple de ce seed.

#     context_kl, context_kr, context_pl, context_pr: pour récupérer le bon voisinage
#     nb_ex_test: nombre d'exemple d'apprentissage a selectionner d'un alignement multiple
#     return:
#         list_example: liste de tuple. Chaque tuple est associé à un exemple qui sera utilisé pour la calcul du score de Brier.
#     """
#     list_pair_seq_name = tuple(dico_seq.keys())
#     # récupération du poid de chaque paire de séquence de dico_seq
#     weights_list = [] 
#     for key in dico_seq:
#         weights_list.append(dico_seq[key]["nbValidPosition"]) 

#     list_pair_name = random.choices(list_pair_seq_name, weights = weights_list, k = nb_ex_test)   # selection des nb_ex_test exemples
#     #print("list_pair_name:\n", list_pair_name)
#     #print("len_list_pair_name:\n", len(list_pair_name))
#     for pair_name in list_pair_name:
#         list_valid_position = dico_seq[pair_name]["validPosition"]  # ensemble des positions valides sans contexte

#         if list_valid_position:   # s'il existe au moins une position valide (a priori)
#             # récupérer la paire dont on sait qu'on va en piocher un exemple (quitte à ce qu'il soit à valeur manquante)
#             pair_seq = dico_seq[pair_name]["SeqAB"]

#             # orienter la paire en 1 couple aléatoire
#             proba = random.uniform(0, 1)                               
#             if proba <= 0.5:
#                 seq_1, seq_2 = pair_seq
#             else:
#                 seq_2, seq_1 = pair_seq   
            
#             # random selection d'une position valide selon le context
#             position_selected = random.choice(list_valid_position)
            
#             len_align = len(seq_1)  # OSD, éventuellement à enlever
 
#             # construction de l'exemple par selection et completion des positions manquantes avec *

#             voisinage_kl = example_shape(seq_1, seq_2, len_align, position_selected, 
#                                          "k", "l", context_kl)
#             voisinage_kr = example_shape(seq_1, seq_2, len_align, position_selected, 
#                                          "k", "r", context_kr)
#             voisinage_pl = example_shape(seq_1, seq_2, len_align, position_selected, 
#                                          "p", "l", context_pl)
#             voisinage_pr = example_shape(seq_1, seq_2, len_align, position_selected, 
#                                          "p", "r", context_pr)
            
#             example_selected = [seq_1[position_selected],  # aa_1
#                                 seq_2[position_selected],  # aa_2
#                                 voisinage_kl,
#                                 voisinage_kr,
#                                 voisinage_pl,
#                                 voisinage_pr]

#             #print("example_selected :", example_selected)


#             #print("position_selected:",position_selected)
#             #print(example_selected)
#             list_example.append(example_selected)
#             #example_selected_count += 1
#             #print("example_selected_count:", example_selected_count)

#     return list_example



# def example_shape(seq_1, seq_2, len_align, position_selected, 
#                   k_or_p, l_or_r, context):
#     """
#     Mise en forme de l'exemple selectionné
#     avec complétion par des * en cas de débordement des indices de l'alignement
#     """
#     voisinage = []
#     index = position_selected
#     # print("seq_1:", seq_1)
#     # print("seq_2:", seq_2)
#     # print("position_selected:", position_selected)

#     if context != 0:
#         if k_or_p == "k":
#             seq = seq_1
#         if k_or_p == "p":
#             seq = seq_2

#         if l_or_r == "l":
#             index -= 1
#             while index >= position_selected - context and index in range(len_align):
#                 voisinage.append(seq[index])
#                 index -= 1

#         else:
#             index += 1
#             while index <= position_selected + context and index in range(len_align):
#                 voisinage.append(seq[index])
#                 index += 1

#         out_of_range = context - len(voisinage)
#         if out_of_range != 0:
#             for i in range(out_of_range):
#                 voisinage.append("*")

#     voisinage = "".join(voisinage)
#     #print("voisinage :", voisinage)
#     return voisinage

    
                

# def multi_random_example_selection(path_folder_seed, path_dico_exemple_complet, path_dico_seq,   # pas encore corrigé
    #                                context_kl, context_kr, context_pl, context_pr):
    # """
    # path_folder_seed: à relire tous les seed max une fois (à voir comment articuler)
    # path_dico_exemple_complet: loader une fois le dico de l'info sur les seed et nb d'exemples à prendre de chacun
    # path_dico_seq: path du folder contenant les dico_seq de chaque seed

    # context_kl: considérer les acides aminés à gauche de l'acide aminé connu jusqu'à context_lk positions à gauche (kl = known left)
    # context_kr: considérer les acides aminés à gauche de l'acide aminé connu jusqu'à context_rk positions à droite (kr = known right)

    # context_pl: considérer les acides aminés à gauche de l'acide aminé à prédire jusqu'à context_lp positions à gauche (pl = to predict left)
    # context_pr: considérer les acides aminés à gauche de l'acide aminé à prédire jusqu'à context_rp positions à droite (pr = to predict right)
    # """
    
    # # initialisation de la liste d'exemples
    # list_example = []

    # # load de dico_exemple une seule fois
    # dico_exemple = np.load(path_dico_exemple_complet, allow_pickle='TRUE').item()

    # # liste des PosixPath des alignements de test
    # files = [x for x in Path(path_folder_seed).iterdir()]
    # nb_files = len(files)
    # #print(nb_files)
    # print("")
    # start = time.time()
    # for file_counter in tqdm(range(nb_files), desc='Selection des ex test /seed',
    #                                ncols= 100, mininterval=60):
    #     file = files[file_counter]
    #     accession_num =  folder.get_accession_number(file)
    #     if accession_num in dico_exemple: 
    #         nb_ex_test = dico_exemple[accession_num]["nbExample"]
    #         if nb_ex_test != 0:

    #             # en loader le dico_seq associé grace à l'accession_num
    #             dico_seq = np.load(f"{path_dico_seq}/{accession_num}.seq.npy", allow_pickle='TRUE').item()
                    
    #             # selection des nb_ex_test exemples dans le seed
    #             list_example = random_example_selection(list_example, dico_seq, 
    #                                                     context_kl, context_kr, context_pl, context_pr, 
    #                                                     nb_ex_test)  

    # print("nb ex test selectionnés :",'{:_}'.format(len(list_example)))

    # end = time.time()
    # diff = end-start
    # items_per_second = nb_files/diff
    # print(f'Selection des ex test: time: {diff:.2f} s | {items_per_second:.2f} it/s')

    # return list_example
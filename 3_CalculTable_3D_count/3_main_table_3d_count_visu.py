################################################################################
#                                  Importations                                #    
################################################################################
import os
import statistics
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.autolayout"] = True

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 


################################################################################
#                                  Variables                                   #    
################################################################################

DATA_RESULT_2D         = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
DATA_RESULT_3D_M2      = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_COUNT"


ALPHABET               = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

liste_tranche          = [(0,10), (10,20), (20,30), (30,40), (40,50), 
                          (50,60),(60,70), (70,80), (80,90), (90,100)]
max_position           = 10

reference              = "origine" 
#reference              = "destination"   

pseudo_compte          = 1

evaluation             = "min_count"      
#evaluation             = "total_count"  
#evaluation             = "not_evaluated"  
################################################################################
#                                  Fonction                                   #    
################################################################################



def descriptionTripletCount(table3d_count, ALPHABET, pseudo_compte):
    """
    table3d_count: dico de la table 3d de comptage a analyser

    renvoie le nombre de triplets comptés au total
            le nombre de triplets non évaluée
            la médiane du nombre d'évaluations des triplets
    """
    # conversion de dico à liste
    list_count_triplet = []

    for aa_1 in ALPHABET:
        for aa_2 in ALPHABET:
            for aa_c in ALPHABET:
                list_count_triplet.append(table3d_count[aa_1][aa_2][aa_c])
    # initialisation du nombre de triplets comptés
    nbre_triplets = 0
    # initialisation du nombre de triplets non évalués
    nbre_triplets_not_evaluated = 0

    for elem in list_count_triplet:
        nbre_triplets += elem
        if elem == pseudo_compte:
            nbre_triplets_not_evaluated += 1 # non évalué ssi l'un des 3 caractères du triplet
                                             # n'est pas dans l'ALPHABET

    min_triplet_evaluation = min(list_count_triplet)
    percentage_not_evaluated = 100*nbre_triplets_not_evaluated/(len(ALPHABET)**3)

    return nbre_triplets, percentage_not_evaluated , min_triplet_evaluation



def descriptionCoupleCount(table2d_count, ALPHABET, pseudo_compte):
    """
    table2d_count: dico de la table 2d de comptage a analyser

    renvoie le nombre de couples comptés au total
            le nombre de couples non évalués
            la médiane du nombre d'évaluations des couples
    """
    # conversion de dico à liste
    list_count_couple= []

    for aa_1 in ALPHABET:
        for aa_2 in ALPHABET:
            list_count_couple.append(table2d_count[aa_1][aa_2])
    # initialisation du nombre de triplets comptés
    nbre_couple = 0
    # initialisation du nombre de triplets non évalués
    nbre_couple_not_evaluated = 0

    for elem in list_count_couple:
        nbre_couple  += elem
        if elem == pseudo_compte:
            nbre_couple_not_evaluated += 1

    min_couple_evaluation = min(list_count_couple)
    percentage_not_evaluated = 100*nbre_couple_not_evaluated/(len(ALPHABET)**2)
    return nbre_couple, percentage_not_evaluated, min_couple_evaluation




plt.figure(figsize=[10,9])
len_list_tranche = len(liste_tranche)

list_position = [i for i in range(-max_position,max_position+1)]

for tranche in liste_tranche:  # tranche de pid
    pid_inf, pid_sup = tranche
    list_total_count = []     # on parle de triplets et de couples ici (tables2d/3d de comptage)
    list_nb_not_evaluated = []
    list_min_count = []

    # chemin des tables 3d count par tranche de pid
    path_table_3d    = f"{DATA_RESULT_3D_M2}/table_3d_count_{pid_inf}_{pid_sup}_{pseudo_compte}"
    
    # importation des tables de comptage 3D de la tranche en position négative
    list_position_neg = [i for i in range(1,max_position+1)][::-1]
    for position_voisin in list_position_neg:
        table3d_comptage = np.load(f"{path_table_3d}/table_3d_count_({-position_voisin},{reference}).npy", allow_pickle='TRUE').item()
        nbre_triplets, percentage_not_evaluated, min_triplet_evaluation = descriptionTripletCount(table3d_comptage, ALPHABET, pseudo_compte)
        list_total_count.append(nbre_triplets)
        list_nb_not_evaluated.append(percentage_not_evaluated)
        list_min_count.append(min_triplet_evaluation)

    # importation de la table de comptage 2D de la tranche
    path_table_2d    = f"{DATA_RESULT_2D}/table_2d_{pid_inf}_{pid_sup}_{pseudo_compte}"

    table2d_comptage = np.load(f"{path_table_2d}/table_2d_count.npy", allow_pickle='TRUE').item()  # à voir comment nommer ces tables2d_count
    nbre_couple, percentage_not_evaluated, min_couple_evaluation = descriptionCoupleCount(table2d_comptage, ALPHABET, pseudo_compte)
    list_total_count.append(nbre_couple)
    list_nb_not_evaluated.append(percentage_not_evaluated)
    list_min_count.append(min_couple_evaluation)

    # importation des tables de comptage 3D de la tranche en position positive
    for position_voisin in range(1,max_position+1):
        table3d_comptage = np.load(f"{path_table_3d}/table_3d_count_({position_voisin},{reference}).npy", allow_pickle='TRUE').item()
        nbre_triplets, percentage_not_evaluated, min_triplet_evaluation = descriptionTripletCount(table3d_comptage, ALPHABET, pseudo_compte)
        list_total_count.append(nbre_triplets)
        list_nb_not_evaluated.append(percentage_not_evaluated)
        list_min_count.append(min_triplet_evaluation)



    if evaluation == "not_evaluated":      #total_count / not_evaluated / median_count
        plt.plot(list_position, list_nb_not_evaluated, label = f"{pid_inf}, {pid_sup}")   # voir après si mettre sur un meme graphe origine/destination
        plt.xticks(range(-max_position,max_position+1))


    if evaluation == "total_count":      #total_count / not_evaluated / median_count
        plt.plot(list_position, list_total_count, label = f"{pid_inf}, {pid_sup}")   # voir après si mettre sur un meme graphe origine/destination
        plt.xticks(range(-max_position,max_position+1))



    if evaluation == "min_count":      #total_count / not_evaluated / median_count
        plt.plot(list_position, list_min_count, label = f"{pid_inf}, {pid_sup}")   # voir après si mettre sur un meme graphe origine/destination
        plt.xticks(range(-max_position,max_position+1))




if evaluation == "not_evaluated": 
    f_str = f"""Position de l'aa contextuel par rapport à l'aa {"d'origine" if reference=="origine" else "de destination"}"""
    plt.xlabel(f_str, fontsize=13)

    plt.ylabel(f"Pourcentage de paramètres non-estimés", fontsize=13)  
    myTitle = f"Evaluation du pourcentage de paramètres non estimés \ndans les tables de comptage 2d et 3d"
    plt.yticks(np.arange(0,60,5))   # ajuste selon les résultats pour améliorer l'affichage
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)                                                    
    plt.title(myTitle, loc='center',  fontsize=20)   # wrap=True
    plt.legend()
    plt.savefig(f"{evaluation}_{reference}.png")


if evaluation == "total_count": 
    f_str = f"""Position de l'aa contextuel par rapport à l'aa {"d'origine" if reference=="origine" else "de destination"}"""
    plt.xlabel(f_str, fontsize=13)

    plt.ylabel(f"Nombre d'exemples d'apprentissage", fontsize=13)
    myTitle = f"Evaluation du nombre d'exemples d'apprentissage \npour le calcul des tables de comptage 2d et 3d"

    plt.yscale('log') 
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)                                                    
    plt.title(myTitle, loc='center', fontsize=20)
    plt.legend()
    plt.savefig(f"{evaluation}_{reference}.png")  # , bbox_inches='tight'     (réduit les bords)


if evaluation == "min_count":     
    f_str = f"""Position de l'aa contextuel par rapport à l'aa {"d'origine" if reference=="origine" else "de destination"}"""
    plt.xlabel(f_str, fontsize=13)

    plt.ylabel(f"Nombre min d'exemples d'apprentissage par paramètre", fontsize=13)   # ou mettre le % car en 0 c'est 400, et ailleurs 8000
    myTitle = f"Evaluation du minimum du nombre d'exemples d'apprentissage par \nparamètre pour le calcul des tables de comptage 2d et 3d"
    plt.yscale('log')  
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)                                                    
    plt.title(myTitle, loc='center',  fontsize=20)
    plt.legend()
    plt.savefig(f"{evaluation}_{reference}.png")

print("visu done")
################################################################################
#                                  Importations                                #    
################################################################################

import matplotlib.pyplot as plt
import os
import argparse
import numpy as np

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 

import brierNeighbour.selection_example as selection_example
import brierNeighbour.brier as brier
import utils.folder as folder



################################################################################
#                                  Variables                                   #    
################################################################################

# variables obligatoires à renseigner lorsqu'on run le script
parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur",
                    type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur",
                    type=int)
parser.add_argument("methode", help="M1 ou M2",
                    type=str)                 
parser.add_argument("nb_voisin", help="_1 si on veut regarder strictement 1 voisin _plus si on veut regarder les voisins jusqu'à une certaine position",
                    type=str)
parser.add_argument("reference", help="k si on prend l'origine comme référence, p si on prend la destination comme référence",
                    type=str)
args = parser.parse_args()


# variables obligatoires à modifier directement dans le script

DATA_MNHN              = f"{file.parents[2]}/MNHN_RESULT/1_DATA"
DATA_TABLE_2D          = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"  # chemin du dossier des tables 2d
DATA_TABLE_3D_PROBA    = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_PROBA_{args.methode}"  # chemin du dossier des tables 3d proba selon la méthode 1
DATA_EXEMPLE_TEST      = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST" # chemin du dossier des sur le pré_traitement des exemples test 
DATA_RESULT            = f"{file.parents[2]}/MNHN_RESULT/5_SCORE_{args.methode}{args.nb_voisin}"  # chemin du dossier des sorties de calcul de score
  
list_standard_aa       = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

nom_folder_fasta_test = "Pfam_split/Pfam_test"

max_position           = 3
nb_test_par_context    = 3
nb_exemple_par_test    = 150 





        
print("#######################################################################")
print("Evaluation de la capacité prédictive d'un acide aminé voisin (méthode1)")
print("#######################################################################")
print("")

################################################################################
#                      Préparation de l'expérience                             #    
################################################################################
  
# origine
if args.reference == "k":
    list_kl          = [(i, 0, 0, 0) for i in range(1,max_position+1)[::-1]] 
    list_non_context = [(0, 0, 0, 0)]
    list_kr          = [(0, i, 0, 0) for i in range(1,max_position + 1)]

    list_context = list_kl  + list_non_context + list_kr 


# destination
if args.reference == "p":
    list_pl          = [(0, 0, i, 0) for i in range(1,max_position+1)[::-1]] 
    list_non_context = [(0, 0, 0, 0)]
    list_pr          = [(0, 0, 0, i) for i in range(1,max_position + 1)]

    list_context = list_pl  + list_non_context + list_pr 
    

list_nb_example = [nb_exemple_par_test]*nb_test_par_context



# création du dossier pour les sorties de score de Brier
# que s'il n'existe pas déjà car on veut faire plusieurs tests sur une meme tranche de pid 
new_file_dico = f"Experience_Brier_{args.pid_inf}_{args.pid_sup}"
path_res_folder = f"{DATA_RESULT}/{new_file_dico}"
isExist = os.path.exists(path_res_folder)
if not isExist:
    os.makedirs(path_res_folder)


# dossier des résultats de score de Brier à référence fixée
path_res_folder_reference = folder.creat_folder(f"{path_res_folder}/Reference_{args.reference}")


# chemin vers les alignements multiples
path_folder_seed = f"{DATA_MNHN}/{nom_folder_fasta_test}"  


# chemin de la table 2d proba de la tranche de pid
path_table_2d_proba = f"{DATA_TABLE_2D}/table_2d_{args.pid_inf}_{args.pid_sup}/table_2d_proba.npy"
# chemin des tables 3d proba de la tranche de pid
path_folder_table_3d_proba_m1    = f"{DATA_TABLE_3D_PROBA}/table_3d_proba_{args.methode}_{args.pid_inf}_{args.pid_sup}"

# chemin des dico d'info sur les seed et seq de pré-traitement des exemples
path_folder_dico_seq = f"{DATA_EXEMPLE_TEST}/Exemple_{args.pid_inf}_{args.pid_sup}/seq"   # chemin ou on veut enregistrer les dico_seq
path_file_dico_seed = f"{DATA_EXEMPLE_TEST}/Exemple_{args.pid_inf}_{args.pid_sup}/seed" # chemin ou on veut enregistrer le dico_seed
path_dico_seed_normalised = f"{path_file_dico_seed}/seed_normalised.npy"  # fraction des ex test a prendre par seed
path_dico_exemple = path_file_dico_seed    
path_dico_exemple_complet = f"{path_file_dico_seed}/exemple_seed.npy"



    




################################################################################
#                                 Expérience                                   #    
################################################################################

list_score_brier  = []
list_position     = []


for context in list_context:
    print("")
    print("##########################")
    print("context :", context)
    print("##########################")

    context_kl, context_kr, context_pl, context_pr = context

    # dans les expériences réalisées avec ce script,
    #  soit les 4 valeurs de context sont nulles, soit une seule est non nulle
    # et on n'étudiera que vers l'origine ou vers la desination (ne pas mélanger les 2)
    # à contraindre ca? ou je serais la seule utilisatrice et donc je sais quels contextes je peux tester?


    total_list_unit_score_brier = []
    for nb_exemple_test in list_nb_example:

        selection_example.example_number_per_seed(path_dico_seed_normalised, nb_exemple_test, path_dico_exemple)
        print("nb ex test demandés:", '{:_}'.format(nb_exemple_test))
        list_example = selection_example.multi_random_example_selection(path_folder_seed, path_dico_exemple_complet, 
                                                                        path_folder_dico_seq, 
                                                                        context_kl, context_kr, context_pl, context_pr)

        #print(list_example)       
        method = args.methode
        if method == "M1":
            if args.nb_voisin == "_1":
                score_brier , list_unit_score_brier = brier.brier_score_m1_1(list_example, 
                                                                     context_kl, context_kr, context_pl, context_pr, 
                                                                     list_standard_aa,
                                                                     path_folder_table_3d_proba_m1, path_table_2d_proba,
                                                                     method)


        if method == "M2":
            if args.nb_voisin == "_1":
                score_brier , list_unit_score_brier = brier.brier_score_m2_1(list_example, 
                                                                     context_kl, context_kr, context_pl, context_pr, 
                                                                     list_standard_aa,
                                                                     path_folder_table_3d_proba_m1, path_table_2d_proba,
                                                                     method)

        for elem in list_unit_score_brier:
            total_list_unit_score_brier.append(elem)


        # origine
        if args.reference == "k" and context != (0,0,0,0):
            if context_kl != 0:
                list_position.append(-context_kl)
                list_score_brier.append(score_brier)
            if context_kr != 0:
                list_position.append(context_kr)
                list_score_brier.append(score_brier)
        
        # destination
        if args.reference == "p" and context != (0,0,0,0):
            if context_pl != 0:
                list_position.append(-context_pl)
                list_score_brier.append(score_brier)
            if context_pr != 0:
                list_position.append(context_pr)
                list_score_brier.append(score_brier)

        # cas non contextuel
        if context == (0,0,0,0):
            list_position.append(0)
            list_score_brier.append(score_brier)


    brier.writeList(total_list_unit_score_brier, f"{path_res_folder_reference}/unit_brier_{context}_{args.reference}")


# nb_test = len(list_nb_example)

# enregistrement des données des graphes
brier.writeList(list_position, f"{path_res_folder_reference}/position_{args.reference}")
brier.writeList(list_score_brier, f"{path_res_folder_reference}/brier_{args.reference}")
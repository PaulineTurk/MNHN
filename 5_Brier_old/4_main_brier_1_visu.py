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


# variables obligatoires à modifier directement dans le script
DATA_RESULT_GLOBAL     = f"{file.parents[2]}/MNHN_RESULT"  # chemin du dossier des sorties 
DATA                   = f"{DATA_RESULT_GLOBAL}/1_DATA" # chemin du dossier des données 
DATA_RESULT            = f"{DATA_RESULT_GLOBAL}/5_SCORE_M1_1"  # chemin du dossier des sorties de calcul de score

nom_folder_fasta_test  = "Pfam_split/Pfam_test"
max_position           = 3
reference              = "p"   # choisir entre k pour "origine" et p pour "destination"
nb_exemple_par_test    = 150
nb_test_par_context    = 3


# liste_tranche = [(0, 10),
#                  (10, 20), 
#                  (20, 30), 
#                  (30, 40),
#                  (40, 50), 
#                  (50, 60), 
#                  (60, 70), 
#                  (70, 80), 
#                  (80, 90),
#                  (90, 100)]

# liste_tranche = [(0, 10)]
# color = "tab:blue"
# liste_tranche = [(10, 20)]
# color = "tab:orange"
# liste_tranche = [(20, 30)]
# color = "tab:green"
# liste_tranche = [(30, 40)]
# color = "tab:red"
# liste_tranche = [(40, 50)]
# color = "tab:purple"
# liste_tranche = [(50, 60)]
# color = "tab:brown"
# liste_tranche = [(60, 70)]
# color = "tab:pink"
#liste_tranche = [(70, 80)]
#color = "tab:grey"
#liste_tranche = [(80, 90)]
#color = "tab:olive"
liste_tranche = [(90, 100)]
color = "tab:cyan"


print("#######################################################################")
print("      Tracé de graphes sur toutes les tranches choisies                ")
print("#######################################################################")
print("")
# création du dossier pour enregistrer les graphes en sortie
path_image = f"{DATA_RESULT}/Graphe_{reference}" 
isExist = os.path.exists(path_image)
if not isExist:
    os.makedirs(path_image)

list_nb_example = [nb_exemple_par_test]*nb_test_par_context
path_folder_seed = f"{DATA}/{nom_folder_fasta_test}"





plt.figure(figsize=[10,9])
len_list_tranche = len(liste_tranche)


for i, tranche in enumerate(liste_tranche):
    pid_inf, pid_sup = tranche
    path_data_graph = f"{DATA_RESULT}/Experience_Brier_{pid_inf}_{pid_sup}/Reference_{reference}"
    list_brier_score = brier.readList(f"{path_data_graph}/brier_{reference}")
    list_position = brier.readList(f"{path_data_graph}/position_{reference}")

    nb_position = len(list_position)
    list_position_err, list_brier_score_mean, list_brier_score_std = brier.estimationError(list_position, list_brier_score, 
                                                                                       nb_position, nb_test_par_context)
    plt.plot(list_position_err, list_brier_score_mean, color=color)
    plt.fill_between(list_position_err, list_brier_score_mean - list_brier_score_std, list_brier_score_mean + list_brier_score_std,
                     alpha=0.5, label = f"{pid_inf}, {pid_sup}", color=color)  #, edgecolor='#CC4F1B', facecolor='#FF9848')


    plt.xticks(range(-max_position,max_position+1))
    #plt.yticks(np.arange(0,1.1,0.1))
    if reference == "k":
        plt.xlabel(f"Position de l'aa contextuel par rapport à l'aa d'origine", fontsize=13)
    if reference == "p":
        plt.xlabel(f"Position de l'aa contextuel par rapport à l'aa de destination", fontsize=13)
    plt.ylabel('Score de Brier (methode 1.1)', fontsize=13)
    myTitle = f"Score de Brier calculé sur {os.path.basename(path_folder_seed)} avec {nb_test_par_context} tests par contexte \net {'{:,}'.format(nb_exemple_par_test).replace(',', ' ')} exemples par test "

    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)                                                    
    plt.title(myTitle, loc='center', wrap=True, fontsize=20)
    plt.legend()
    
    title = f"Score de Brier calculé sur {os.path.basename(path_folder_seed)}"
    if reference == "k":
        plt.savefig(f"{path_image}/{title}_origine_{pid_inf}_{pid_sup}.png")
    if reference == "p":
        plt.savefig(f"{path_image}/{title}_destination_{pid_inf}_{pid_sup}.png")
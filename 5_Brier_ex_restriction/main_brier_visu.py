"""
Visualisation of the Brier Score tests
"""

################################################################################
#                                  Importations                                #
################################################################################

import os
import sys
from pathlib import Path
# import argparse
# import numpy as np
import matplotlib.pyplot as plt



# import brierNeighbour.selection_example as selection_example
import brierNeighbour.brier as brier
# import utils.folder as folder

file = Path(__file__).resolve()
sys.path.append(file.parents[0])





################################################################################
#                                  Variables                                   #
################################################################################
METHODE = "uni"

# variables obligatoires à modifier directement dans le script

# chemin du dossier des sorties
DATA_RESULT_GLOBAL = f"{file.parents[2]}/MNHN_RESULT"

# chemin du dossier des données
DATA = f"{DATA_RESULT_GLOBAL}/1_DATA"

# chemin du dossier des sorties de calcul de score
DATA_RESULT = f"{DATA_RESULT_GLOBAL}/5_SCORE_{METHODE}_ex_restriction"


NOM_FOLDER_FASTA_TEST = "Pfam_split/Pfam_test"
MAX_POSITION = 10
# REFERENCE = "origine"
REFERENCE = "destination"
NB_EXEMPLE_PAR_TEST = 1_500_000
NB_TEST_PAR_CONTEXT = 5

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

# liste_tranche = [(0, 10),
#                  (70, 80),
#                  (80, 90),
#                  (90, 100)]

liste_tranche = [(0, 10)]
COLOR = "tab:blue"
# liste_tranche = [(10, 20)]
# COLOR = "tab:orange"
# liste_tranche = [(20, 30)]
# COLOR = "tab:green"
# liste_tranche = [(30, 40)]
# COLOR = "tab:red"
# liste_tranche = [(40, 50)]
# COLOR = "tab:purple"
# liste_tranche = [(50, 60)]
# COLOR = "tab:brown"
# liste_tranche = [(60, 70)]
# COLOR = "tab:pink"
# liste_tranche = [(70, 80)]
# COLOR = "tab:grey"
# liste_tranche = [(80, 90)]
# COLOR = "tab:olive"
# liste_tranche = [(90, 100)]
# COLOR = "tab:cyan"


print("#######################################################################")
print("      Tracé de graphes sur toutes les tranches choisies                ")
print("#######################################################################")
print("")
# création du dossier pour enregistrer les graphes en sortie
path_image = f"{DATA_RESULT}/Graphe_{REFERENCE}"
IS_EXIST = os.path.exists(path_image)
if not IS_EXIST:
    os.makedirs(path_image)

list_nb_example = [NB_EXEMPLE_PAR_TEST]*NB_TEST_PAR_CONTEXT
path_folder_seed = f"{DATA}/{NOM_FOLDER_FASTA_TEST}"





plt.figure(figsize=[10,9])

for i, tranche in enumerate(liste_tranche):
    pid_inf, pid_sup = tranche
    path_data_graph = f"{DATA_RESULT}/Experience_Brier_{pid_inf}_{pid_sup}/Reference_{REFERENCE}"
    list_brier_score = brier.readList(f"{path_data_graph}/brier_{REFERENCE}")
    list_position = brier.readList(f"{path_data_graph}/position_{REFERENCE}")

    nb_position = len(list_position)
    list_position_err, list_brier_score_mean, list_brier_score_std = brier.estimationError(
                         list_position, list_brier_score, nb_position, NB_TEST_PAR_CONTEXT)

    plt.plot(list_position_err, list_brier_score_mean, color=COLOR)
    plt.fill_between(list_position_err, list_brier_score_mean - list_brier_score_std,
                     list_brier_score_mean + list_brier_score_std,
                     alpha=0.5, label = f"{pid_inf}, {pid_sup}", color=COLOR)

    plt.xticks(range(-MAX_POSITION,MAX_POSITION+1))
    #plt.yticks(np.arange(0,1.1,0.1))
    f_str = f"""Position de l'aa contextuel par rapport à l'aa \
    {"d'origine" if REFERENCE=="origine" else "de destination"}"""
    plt.xlabel(f_str, fontsize=13)
    plt.ylabel(f'Score de Brier {METHODE}', fontsize=13)
    NAME_FOLDER_MSA = os.path.basename(path_folder_seed)
    myTitle = f"Score de Brier calculé sur {NAME_FOLDER_MSA} avec {NB_TEST_PAR_CONTEXT} tests par contexte\n\
               {NB_EXEMPLE_PAR_TEST:,} exemples par test"   # à rechanger
    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)
    plt.title(myTitle, loc='center', fontsize=20)
    plt.legend()
    title = f"Score de Brier calculé sur {os.path.basename(path_folder_seed)}"
    plt.savefig(f"{path_image}/{title}_{REFERENCE}_{pid_inf}_{pid_sup}.png")

"""
Preprocessing of data selection
"""
################################################################################
#                                  Importations                                #
################################################################################

import os
import argparse
# import numpy as np
# import matplotlib.pyplot as plt


import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

import data_exemple_test.dataexampletestfonction as dataexampletestfonction
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
args = parser.parse_args()


# variables obligatoires à modifier directement dans le script

# chemin du dossier des données
DATA = f"{file.parents[2]}/MNHN_RESULT/1_DATA"

# chemin du dossier des sorties du pré-traitement de selection des ex test
DATA_RESULT  = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST"
ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

NOM_FOLDER_FASTA_TEST = "Pfam_split/Pfam_test"
NOM_FOLDER_PID = "PID"


print("#######################################################################")
print("     Pré-traitement pour la selection des exemples de test             ")
print("#######################################################################")
print("")

path_folder_seed = f"{DATA}/{NOM_FOLDER_FASTA_TEST}"   # chemin vers les seeds tests
path_folder_pid = f"{DATA}/{NOM_FOLDER_PID}"  # chemin vers les pid de seeds

# création du dossier pour les dico de pré_traitement des exemples test
# que s'il n'existe pas déjà
path_res_folder = f"{DATA_RESULT}"
IS_EXIST = os.path.exists(path_res_folder)
if not IS_EXIST:
    os.makedirs(path_res_folder)


# tout ce qui suit est relatif à une tranche de pid :

# création du dossier pour les dico de pré-traitement
new_folder_dico = f"Exemple_{args.pid_inf}_{args.pid_sup}"
path_res_folder = folder.creat_folder(f"{DATA_RESULT}/{new_folder_dico}")

 # chemin ou on veut enregistrer l'information par alignement multiple
path_folder_dico_seq = f"{path_res_folder}/seq"

 # chemin ou on veut enregistrer l'information sur tous les alignements multiples
path_file_dico_seed = f"{path_res_folder}/seed"

# calcul du dico_seed et des dico_seq
dataexampletestfonction.multi_seeds_selection(path_folder_seed, path_folder_pid,
                                              args.pid_inf, args.pid_sup, ALPHABET,
                                              path_folder_dico_seq, path_file_dico_seed)

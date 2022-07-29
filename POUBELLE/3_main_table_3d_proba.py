################################################################################
#                                  Importations                                #    
################################################################################
import os
import argparse
import numpy as np

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 

import table3d_proba.table3dprobafonction as table3dprobafonction
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
parser.add_argument("max_position", help="max_position voisin",
                    type=int)
parser.add_argument("kp", help="'origine' ou 'destination'",
                    type=str)
args = parser.parse_args()


# variables obligatoires à modifier directement dans le script

DATA_COUNT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_COUNT" # chemin du dossier des tables de comptage
DATA_RESULT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_PROBA"  # chemin du dossier des sorties

ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


pseudo_compte = 1

print("#######################################################################")
print("    Calcul table 3d proba à partir de la table 3d count associée       ")
print("#######################################################################")

kp_SeqChoice = args.kp
delay_num = args.max_position
print(f"\nTable 3d proba ({delay_num},{kp_SeqChoice}) de la tranche de pid [{args.pid_inf},{args.pid_sup}]\n")
                    

# creation du dossier d'enregistrement des tables 3d proba
# s'il n'existe pas déjà
path_table_3d = f"{DATA_RESULT}"
isExist = os.path.exists(path_table_3d)
if not isExist:
    os.makedirs(path_table_3d)


# creation du dossier d'enregistrement des tables 3d proba selon la tranche de pid
# s'il n'existe pas déjà
path_table_3d_par_tranche_de_pid = f"{DATA_RESULT}/table_3d_proba_{args.pid_inf}_{args.pid_sup}_{pseudo_compte}"
isExist = os.path.exists(path_table_3d_par_tranche_de_pid)
if not isExist:
    os.makedirs(path_table_3d_par_tranche_de_pid)


# chargement de la table 3d count
path_table_3d_count = f"{DATA_COUNT}/table_3d_count_{args.pid_inf}_{args.pid_sup}_{pseudo_compte}"
table_3d_count = np.load(f"{path_table_3d_count}/table_3d_count_({delay_num},{kp_SeqChoice}).npy", allow_pickle='TRUE').item()


table_3d_proba_m2, table_3d_proba_m2= table3dprobafonction.table_3d_proba(table_3d_count, delay_num, kp_SeqChoice, 
                                                                          ALPHABET,
                                                                          path_table_3d_par_tranche_de_pid)
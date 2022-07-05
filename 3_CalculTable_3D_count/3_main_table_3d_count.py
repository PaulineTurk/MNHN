################################################################################
#                                  Importations                                #    
################################################################################
import os
import argparse

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 

import table3d_count.table3dcountfonction as table3dcountfonction
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
parser.add_argument("position_relative_voisin", help="position du voisin par rapport la position de référence",
                    type=int)
parser.add_argument("seq_reference", help="choisir la séquence de référence pour le voisin 'origine' ou 'destination'",
                    type=str)
args = parser.parse_args()


# variables obligatoires à modifier directement dans le script
DATA                   = f"{file.parents[2]}/MNHN_RESULT/1_DATA" # chemin du dossier des données 
DATA_RESULT            = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_COUNT"  # chemin du dossier des sorties 

ALPHABET               = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

nom_folder_fasta_train = "Pfam_split/Pfam_train"
nom_folder_PID         = "PID"
pseudo_compte          = 1



################################################################################
#                         Calcul des tables 3d de comptage                     #    
################################################################################
print("#######################################################################")
print("                         Calcul table 3d de comptage                   ")
print("#######################################################################")
path_folder_fasta      = f"{DATA}/{nom_folder_fasta_train}"
path_folder_fasta      = f"{DATA}/{nom_folder_fasta_train}"
path_folder_pid        = f"{DATA}/{nom_folder_PID}"


# création du dossier globale d'enregistrement des tables 3d,
# que s'il n'existe pas
isExist = os.path.exists(DATA_RESULT)
if not isExist:
    os.makedirs(DATA_RESULT)

# création du dossier globale d'enregistrement des tables 3d,
# que s'il n'existe pas
path_table_3d_count = f"{DATA_RESULT}/table_3d_count_{args.pid_inf}_{args.pid_sup}_{pseudo_compte}"
isExist = os.path.exists(path_table_3d_count)
if not isExist:
    os.makedirs(path_table_3d_count)



delay_num = args.position_relative_voisin
RefSeqChoice = args.seq_reference
print(f"\n{delay_num}, {RefSeqChoice}\n")
table3dcountfonction.multi_table_3d_count(path_folder_fasta, path_folder_pid, 
                                          delay_num, RefSeqChoice, 
                                          ALPHABET, args.pid_inf, args.pid_sup,
                                          path_table_3d_count,
                                          pseudo_compte)
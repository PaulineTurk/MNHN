################################################################################
#                                  Importations                                #    
################################################################################
import os
import argparse


import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 

import table2d.table2dfonction as table2dfonction
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
DATA                   = f"{file.parents[2]}/MNHN_RESULT/1_DATA" # chemin du dossier des données 
DATA_RESULT            = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"  # chemin du dossier des sorties 

list_standard_aa       = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

nom_folder_fasta_train = "Pfam_split/Pfam_train"
nom_folder_PID         = "PID"

pseudo_compte          = 1

################################################################################
#                       Compte des acides aminés standards                     #
#                   et des couples d'acides aminés standards                   #    
################################################################################

# création du dossier globale d'enregistrement des tables 2d,
# que s'il n'existe pas
isExist = os.path.exists(DATA_RESULT)
if not isExist:
    os.makedirs(DATA_RESULT)


print("#######################################################################")
print("                             Calcul table 2d                           ")
print("#######################################################################")
path_folder_fasta      = f"{DATA}/{nom_folder_fasta_train}"
path_folder_pid        = f"{DATA}/{nom_folder_PID}"

print("")
print(f"table_2d ({args.pid_inf},{args.pid_sup})")


# création du dossier d'enregistrement des tables 2d pour une tranche de pid,
path_folder_Result = f"{DATA_RESULT}/table_2d_{args.pid_inf}_{args.pid_sup}_{pseudo_compte}"
path_folder_Result = folder.creat_folder(path_folder_Result)


################################################################################
#                     Calcul table 2d de comptage                              #    
################################################################################

count_AA, nb_AA, count_coupleAA, nb_coupleAA = table2dfonction.multi_count_for_table_2d(path_folder_fasta, path_folder_pid, 
                                                                                     list_standard_aa, args.pid_inf, args.pid_sup,
                                                                                     path_folder_Result,
                                                                                     pseudo_compte)


################################################################################
#                     Calcul table 2d de fréquences                            #    
################################################################################
                                                                                 
freq_AA, freq_coupleAA = table2dfonction.freq_for_table_2d(count_AA, nb_AA, count_coupleAA, nb_coupleAA, path_folder_Result)



################################################################################
#              Calcul table 2d de score (formule de BLOSUM)                    #    
################################################################################

table_2d_score = table2dfonction.table_2d_score(freq_AA, freq_coupleAA, path_folder_Result, scale_factor = 2)


################################################################################
#                          Calcul table 2d de proba                            #    
################################################################################

cond_proba = table2dfonction.table_2d_conditional_proba(freq_AA, freq_coupleAA, path_folder_Result)


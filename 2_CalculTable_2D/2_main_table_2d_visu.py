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



path_folder_fasta      = f"{DATA}/{nom_folder_fasta_train}"
path_folder_pid        = f"{DATA}/{nom_folder_PID}"

pseudo_compte          = 1

print("#######################################################################")
print("                       Visualisation table 2d                          ")
print("#######################################################################")


################################################################################
#           Chemin vers les tables 2d de la tranche {pid_inf},{pid_sup}        #    
################################################################################

path_folder_Result_par_tranche_de_pid = f"{DATA_RESULT}/table_2d_{args.pid_inf}_{args.pid_sup}_{pseudo_compte}"  
path_result_graphe = folder.creat_folder(f"{path_folder_Result_par_tranche_de_pid}/graphe")

################################################################################
#                       Visualisation table 2d de score                        #    
################################################################################
print(f"\ntable_2d_score ({args.pid_inf},{args.pid_sup}) :\n")
table_2d_score =  np.load(f"{path_folder_Result_par_tranche_de_pid}/table_2d_score.npy", allow_pickle='TRUE').item()  # a vérifier nom
table2dfonction.table_2d_visualisation_transposition(table_2d_score)
name_folder_fasta =  os.path.basename(path_folder_fasta)
title_heatmap = f"Heatmap de la table_2d de scores [{args.pid_inf},{args.pid_sup}] calculée sur {name_folder_fasta}"
table2dfonction.table_2d_heatmap(table_2d_score, path_result_graphe, title_heatmap, 
                                 size_annot = 5)


################################################################################
#    Visualisation comparaison table 2d de score et BLOSUM_62 de référence     #    
################################################################################
pid_inf_ref = 62
matrix_diff, pid_inf_ref, average_diff = table2dfonction.table_2d_difference(table_2d_score, pid_inf_ref)
title_heatmap = f"Heatmap des différences entre la table_2d de scores [{args.pid_inf},{args.pid_sup}] et la Blosum_{pid_inf_ref} de référence\nLa différence moyenne de score est de : {average_diff}"
table2dfonction.table_2d_heatmap(matrix_diff, path_result_graphe, title_heatmap, 
                                 size_annot = 5)


################################################################################
#                   Visualisation table 2d de proba                            #    
################################################################################
print(f"\ntable_2d_proba ({args.pid_inf},{args.pid_sup}) :\n")
table_2d_proba =  np.load(f"{path_folder_Result_par_tranche_de_pid}/table_2d_proba.npy", allow_pickle='TRUE').item()  # a vérifier nom

table2dfonction.table_2d_visualisation_transposition(table_2d_proba)
table2dfonction.sum_line_transposition(table_2d_proba)
title_heatmap = f"Heatmap de la table_2d de probabilités conditionnelles [{args.pid_inf},{args.pid_sup}] calculée sur {name_folder_fasta}"
table2dfonction.table_2d_heatmap(table_2d_proba, path_result_graphe, title_heatmap,
                                 size_annot = 3)


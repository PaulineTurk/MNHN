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
DATA_RESULT            = f"{DATA_RESULT_GLOBAL}/4_SCORE_M1_1"  # chemin du dossier des sorties de calcul de score

pid_inf, pid_sup = 90, 100
context = (0,0,0,2)
new_file_dico = f"Experience_Brier_{pid_inf}_{pid_sup}"
path_res_folder = f"{DATA_RESULT}/{new_file_dico}"



list_brier_score = brier.readList("/home/pauline/Bureau/MNHN_RESULT/5_SCORE_M1_1/Experience_Brier_90_100_150/Reference_p/unit_brier_(0, 0, 1, 0)_p")
print(len(list_brier_score))
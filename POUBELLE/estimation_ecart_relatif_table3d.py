################################################################################
#                                  Importations                                #    
################################################################################
import os
import statistics
import numpy as np
import matplotlib.pyplot as plt

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[0]) 
import utils.folder as folder

################################################################################
#                                  Variables                                   #    
################################################################################

DATA_RESULT            = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D_COUNT"  # chemin du dossier des sorties 

list_standard_aa       = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


liste_tranche = [(0, 10),
                 (10, 20), 
                 (20, 30), 
                 (30, 40),
                 (40, 50), 
                 (50, 60), 
                 (60, 70), 
                 (70, 80), 
                 (80, 90),
                 (90, 100)]
max_position  = 10

pseudo_compte = 1

#operation = "max"
#operation = "mean"
operation = "median"


################################################################################
#                         Fonctions tables 3d count                            #    
################################################################################


def ecartRelatifkp(path_tables_3d, position, pid_inf, pid_sup, list_standard_aa):
    """
    position: décalage de position, négatif si à gauche, positif si à droite
    """

    table3d_k = np.load(f"{path_tables_3d}/table_3d_count_{pid_inf}_{pid_sup}_{pseudo_compte}/table_3d_count_({position},origine).npy", allow_pickle='TRUE').item()
    table3d_p = np.load(f"{path_tables_3d}/table_3d_count_{pid_inf}_{pid_sup}_{pseudo_compte}/table_3d_count_({position},destination).npy", allow_pickle='TRUE').item()

    # intialisation de la liste des écart-relatif absolu terme à terme table3d_k - table3d_p
    ecart_relatif_k_p = []
    # initialisation du max
    maximum = 0
    # calcul de cette liste
    for aa_1 in list_standard_aa:
        for aa_2 in list_standard_aa:
            for aa_c in list_standard_aa:
                if table3d_k[aa_1][aa_2][aa_c]!= pseudo_compte and table3d_p[aa_1][aa_2][aa_c]!= pseudo_compte: # ne considérer que les paramètres qui ont été estimés
                    ecart_relatif = abs(table3d_k[aa_1][aa_2][aa_c] - table3d_p[aa_1][aa_2][aa_c])/statistics.mean([table3d_k[aa_1][aa_2][aa_c], table3d_p[aa_1][aa_2][aa_c]])
                    ecart_relatif_k_p.append(ecart_relatif)
                    if maximum < ecart_relatif:
                        maximum = ecart_relatif
                        position_max = (aa_1, aa_2, aa_c)
    # print(position, pid_inf, pid_sup)
    # print("triplet_max:", position_max)
    #print("ecart_relatif_k_p: ", ecart_relatif_k_p)
    #print("len ecart_relatif_k_p:", len(ecart_relatif_k_p))
    return ecart_relatif_k_p


def max_ecart_relatif(max_ecart_relatif_par_tranche, list_ecart_relatif):
    max_ecart_relatif_par_tranche.append(max(list_ecart_relatif))
    return max_ecart_relatif_par_tranche

def mean_ecart_relatif(mean_ecart_relatif_par_tranche, list_ecart_relatif):
    mean_ecart_relatif_par_tranche.append(statistics.mean(list_ecart_relatif))
    return mean_ecart_relatif_par_tranche

def median_ecart_relatif(median_ecart_relatif_par_tranche, list_ecart_relatif):
    median_ecart_relatif_par_tranche.append(statistics.median(list_ecart_relatif))
    return median_ecart_relatif_par_tranche

################################################################################
#                        Comparaison des tables 3d count                       #    
################################################################################


dico_fn_comparaison = {"max": max_ecart_relatif, 
                       "mean": mean_ecart_relatif,
                       "median": median_ecart_relatif}




plt.figure(figsize=[10,9])
len_list_tranche = len(liste_tranche)
path_tables_3d = DATA_RESULT
list_position = [i for i in range(-max_position,max_position+1) if i!=0]

for i, tranche in enumerate(liste_tranche):  # tranche de pid
    pid_inf, pid_sup = tranche
    ecart_relatif_par_tranche = []
    for position in list_position:
        list_ecart_relatif = ecartRelatifkp(path_tables_3d, position, pid_inf, pid_sup, list_standard_aa)
        ecart_relatif_par_tranche = dico_fn_comparaison[operation](ecart_relatif_par_tranche, list_ecart_relatif)

    plt.plot(list_position, ecart_relatif_par_tranche, label = f"{pid_inf}, {pid_sup}")
    plt.xticks(range(-max_position,max_position+1))
    plt.xlabel(f"Position de l'aa contextuel par rapport à l'aa de référence", fontsize=13)
    plt.ylabel(f"{operation}(abs(Ecart-relatif))", fontsize=13)
    myTitle = f"Evaluation de la différence entre table3d origine et table3d destination -{operation}"

    plt.grid(color='lightgrey', linestyle='--', linewidth=0.5)                                                    
    plt.title(myTitle, loc='center', wrap=True, fontsize=20)
    plt.legend()
    plt.savefig(f"Difference_origine_destination_table3d_count_{operation}.png")

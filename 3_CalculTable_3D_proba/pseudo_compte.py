"""
Create 3d_tables with different weight of 2d_table
"""
import sys
import argparse
from pathlib import Path
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("pid_inf", help="pourcentage d'identité inférieur",
                    type=int)
parser.add_argument("pid_sup", help="pourcentage d'identité supérieur",
                    type=int)
parser.add_argument("delay_num", help="max_position voisin",
                    type=int)
parser.add_argument("origine_destination", help="'origine' ou 'destination'",
                    type=str)
args = parser.parse_args()

file = Path(__file__).resolve()
sys.path.append(file.parents[1])

PSEUDO_COMPTE = 1
ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
PATH_TABLE_3D_PROBA = ""
DATA_3D_COUNT = f"{file.parents[2]}/MNHN_RESULT/3_TABLE_3D"
DATA_2D_COUNT = f"{file.parents[2]}/MNHN_RESULT/2_TABLE_2D"
WEIGHT = 0.1



# Fonctions
def table_3d_remove_1(table_3d_count, alphabet, value):
    """
    Remove a value from each cell of the 3d_table_count
    """
    for aa_origine in alphabet:
        for aa_destination in alphabet:
            for aa_context in alphabet:
                table_3d_count[aa_origine][aa_destination][aa_context] -= value
    return table_3d_count


def table_3d_proba(table_3d_count, alphabet):
    """
    Compute 3d_table_proba
    """
    intra_couple_count = {}
    for aa_k in alphabet:
        intra_couple_count[aa_k] = {}
        for aa_p in alphabet:
            intra_couple_count[aa_k][aa_p] = 0
            for aa_c in alphabet:
                intra_couple_count[aa_k][aa_p] += table_3d_count[aa_k][aa_p][aa_c]

    table_3d_proba = {}
    for aa_k in alphabet:
        table_3d_proba[aa_k] = {}
        for aa_p in alphabet:
            table_3d_proba[aa_k][aa_p] = {}
            for aa_c in alphabet:
                if intra_couple_count[aa_k][aa_p] != 0:
                    table_3d_proba[aa_k][aa_p][aa_c] = (table_3d_count[aa_k][aa_p][aa_c]) / (intra_couple_count[aa_k][aa_p])  
                else:
                    table_3d_proba[aa_k][aa_p][aa_c] = 0
    return table_3d_proba


# importer table 3D de comptage
PATH_TABLE_3D_COUNT = f"{DATA_3D_COUNT}/table_3d_count_{args.pid_inf}_{args.pid_sup}_{PSEUDO_COMPTE}"
TABLE_3D_COUNT = np.load(f"{PATH_TABLE_3D_COUNT}/table_3d_count_({args.delay_num},{args.origine_destination}).npy",
                           allow_pickle='TRUE').item()

# en soustraire le pseudo_compte si nécessaire
TABLE_3D_COUNT_SS_PSEUDO_COMTPE = table_3d_remove_1(ALPHABET, TABLE_3D_COUNT, PSEUDO_COMPTE)

# en déduire les tables de proba sans pseudo_compte
TABLE_3D_PROBA = table_3d_proba(TABLE_3D_COUNT_SS_PSEUDO_COMTPE, ALPHABET)

# calculer des tables 3d proba selon une contribution alpha de la table_2d_proba
# pas grave si je garde le pseudo_compte à 1 ici (cela assure qu'on n'ajoutera pas de 0 dans les
# tables 3d proba)
def pseudo_weight_modulation(table_3d_proba, table_2d_proba, weight, alphabet):
    """Weight the table_3d_proba with the table_2d_proba

    Args:
        table_3d_proba (dico): table_3d_proba
        table_2d_proba (dico): table_3d_proba
        weight (float): contribution ajout
        alphabet (list): list of the valid characters
    """
    for aa_origine in alphabet:
        for aa_destination in alphabet:
            for aa_context in alphabet:
                prior = table_2d_proba[aa_origine][aa_destination]
                table_3d_proba[aa_origine][aa_destination][aa_context] += weight*prior
    return table_3d_proba

    
TABLE_3D_PROBA_MODULATED = pseudo_weight_modulation(TABLE_3D_PROBA, table_2d_proba, WEIGHT, ALPHABET)
# je n'ai pas renormalisé après cette étape, est-ce nécessaire ????????????????

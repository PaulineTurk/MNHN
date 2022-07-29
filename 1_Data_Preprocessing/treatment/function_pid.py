# IMPORTS

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[1])


# FUNCTIONS

def pid(seq_1: str, seq_2: str,
        alphabet: list, len_align: int):
    """Compute the percentage of identity between 2 sequences

    Args:
        seq_1 (str): first sequence
        seq_2 (str): second sequence
        alphabet (list): list of character included
        len_align: number of characters in seq_1 and seq_2

    Returns:
        float: percentage of identity
    """
    pid = 0
    n_included_character_seq_1 = 0
    n_included_character_seq_2 = 0

    for indice_aa in range(len_align):
        inclusion_check = 0
        if seq_1[indice_aa] in alphabet:
            n_included_character_seq_1 += 1
            inclusion_check += 1
        if seq_2[indice_aa] in alphabet:
            n_included_character_seq_2 += 1
            inclusion_check += 1
            
        if seq_1[indice_aa] == seq_2[indice_aa] and inclusion_check == 2:
            pid += 1

    return 100*pid/min(n_included_character_seq_1, n_included_character_seq_2)

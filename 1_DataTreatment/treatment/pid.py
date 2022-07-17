# IMPORTS

import numpy as np
from tqdm import tqdm
import time

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.fastaReader as fastaReader
import utils.folder as folder



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




def save_pid(path_folder_fasta: str, path_folder_pid: str, alphabet: list):
    """Compute the percentage of identity for each couple of sequences in a
       folder of multi-sequence alignments

    Args:
        path_folder_fasta (str): path od data to compute pid on
        path_folder_pid (str): path where the pid info is stored
        alphabet (list): list of character included
    """

    folder.creat_folder(path_folder_pid)

    files = [x for x in Path(path_folder_fasta).iterdir()]
    n_files = len(files)

    start = time.time()
    
    for file_counter in tqdm(range(n_files), desc = "pid", mininterval=60):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        path_file_pid = f"{path_folder_pid}/{accession_num}.pid"
        liste_seq = fastaReader.read_multi_fasta(file)
        pid_couple = {}
        n_seq = len(liste_seq)
        len_align = int(len(liste_seq[0][1]))
        for i in range(n_seq):
            pid_couple[liste_seq[i][0]] = {}
        
        for i in range(n_seq):
            for j in range(i, n_seq):
                current_pid = pid(liste_seq[i][1], liste_seq[j][1],
                                  alphabet, len_align)
                pid_couple[liste_seq[i][0]][liste_seq[j][0]] = current_pid
                pid_couple[liste_seq[j][0]][liste_seq[i][0]] = current_pid
        
        np.save(path_file_pid, pid_couple)
    
    end = time.time()
    diff = end - start
    print(f"PID: time {diff} s")

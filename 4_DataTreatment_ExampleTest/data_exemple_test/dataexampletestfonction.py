# IMPORTS
import os
import time
import sys
from pathlib import Path
import numpy as np
from tqdm import tqdm
 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])

import utils.fastaReader as fastaReader 
import utils.folder as folder


def one_seed_selection(accession_num, seed,
                       pid, pid_inf, pid_sup,
                       n_non_info_seed,  n_valid_aa_couple_global,
                       dico_seed,
                       alphabet):

    len_align = len(seed[0][1])
    n_seq = len(seed)

    dico_seq = {}
    n_valid_aa_couple_seed = 0

    for i in range(n_seq):
        name_A, seq_A = seed[i]
        for j in range(i+1, n_seq):
            name_B, seq_B = seed[j]

            if pid_inf <= pid[name_A][name_B] < pid_sup:
                list_valid_index = []
                for index in range(len_align):
                    if all(x in alphabet for x in [seq_A[index],seq_B[index]]):
                        list_valid_index.append(index)
                n_valid_aa_couple = len(list_valid_index)
                        
                if n_valid_aa_couple != 0:
                    dico_seq[(name_A, name_B)] = {}
                    dico_seq[(name_A, name_B)]["N_VALID_AA_COUPLE"] = n_valid_aa_couple
                    dico_seq[(name_A, name_B)]["SEQA_SEQB"] = (seq_A, seq_B)
                    dico_seq[(name_A, name_B)]["LIST_VALID_INDEX"] = list_valid_index
                                          
                n_valid_aa_couple_seed += n_valid_aa_couple

    if n_valid_aa_couple_seed == 0: # not informative seed
        n_non_info_seed += 1

    else:
        dico_seed[accession_num] = {}
        dico_seed[accession_num]["N_VALID_AA_COUPLE"] = n_valid_aa_couple_seed
        dico_seed[accession_num]["LEN_ALIGN"] = len_align

    n_valid_aa_couple_global += n_valid_aa_couple_seed

    return dico_seed, n_non_info_seed, n_valid_aa_couple_global, dico_seq


def multi_seeds_selection(path_folder_seed, path_folder_pid,
                          pid_inf, pid_sup, alphabet,
                          path_folder_dico_seq, path_folder_dico_seed):
    # """
    # path_folder_dico_seq: dossier créé s'il n'existe pas, effacé puis recréé sinon.
    #                       Les dictionnaires de seed y sont enregistrés.
    # path_file_dico_seed: dossier créé s'il n'existe pas, effacé puis recréé sinon.
    #                      Le dictionnaire sur les seeds y est enregistré
    #                      ainsi que sa version normalisée par la somme de tous les comptes.
    # """

    # GLOBAL COUNTERS INITIALISATION
    dico_seed = {}
    n_non_info_seed = 0
    n_valid_aa_couple_global = 0

    files = [x for x in Path(path_folder_seed).iterdir()]
    n_files = len(files)
    print("N_MSA: ", n_files)
    start = time.time()
    for file_counter in tqdm(range(n_files), desc='VALID MSA SELECTION',
                                   ncols= 100, mininterval=60):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        pid = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()
        seed = fastaReader.read_multi_fasta(file)

        dico_seed, n_non_info_seed, n_valid_aa_couple_global, dico_seq = one_seed_selection(
                                                                            accession_num, seed,
                                                                            pid, pid_inf, pid_sup,
                                                                            n_non_info_seed,
                                                                            n_valid_aa_couple_global,
                                                                            dico_seed,
                                                                            alphabet)
        
        # DICO SEQ SAVE
        path_file_dico_seq = f"{path_folder_dico_seq}/{accession_num}.seq"
        np.save(path_file_dico_seq, dico_seq)

    # DICO SEED SAVE
    path_save_file_dico_seed = f"{path_folder_dico_seed}/seed"
    np.save(path_save_file_dico_seed, dico_seed)

    # DICO SEED NORMALISATION
    dico_seed_normal = {}
    for key in dico_seed:
        dico_seed_normal[key] = {}
        dico_seed_normal[key]["WEIGHT_SEED"] = dico_seed[key]["N_VALID_AA_COUPLE"]/n_valid_aa_couple_global
        dico_seed_normal[key]["LEN_ALIGN"] = dico_seed[key]["LEN_ALIGN"]
        
    # DICO SEED NORMALISED SAVE
    path_save_file_dico_seed_normalised = f"{path_folder_dico_seed}/seed_normal"
    np.save(path_save_file_dico_seed_normalised, dico_seed_normal)

    n_valid_seed = len(dico_seed)
    print(f"\nN_VALID, N_INVALID SEED: {n_valid_seed}, {n_non_info_seed}")
    end = time.time()
    diff = end - start
    items_per_second = n_files/diff
    print(f'TIME PRE-PROCESSING TEST EXAMPLES: {diff:.2f}s | {items_per_second:.2f}it/s')

    return n_valid_seed, n_non_info_seed, n_valid_aa_couple_global

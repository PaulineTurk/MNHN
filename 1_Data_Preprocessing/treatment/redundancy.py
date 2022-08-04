# IMPORTS

import sys
import time
from pathlib import Path
import numpy as np
from tqdm import tqdm

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.fastaReader as fastaReader
import utils.folder as folder


# FUNCTIONS

def len_seq_corrected(seq, alphabet):
    """
    Return the number of residus in seq that are included in alphabet
    """
    len_seq_corrected = 0
    for aa in seq:
        if aa in alphabet:
            len_seq_corrected += 1
    return len_seq_corrected



def clustering_non_redundant(liste_seq, file_seq_non_redondant,
                             alphabet, pid_cluster, data_pid):
    """
    Return a partition of liste_seq of sequences with a percentage of identity greater or equal than pid_cluster
    """
    cluster = {}
    if liste_seq:
        name_0, seq_0 = liste_seq[0] 
        len_seq_real_0 = len_seq_corrected(seq_0, alphabet)
        cluster[0] = [(name_0, seq_0, len_seq_real_0)]

        for name_1, seq_1 in liste_seq:
            len_seq_real_1 = len_seq_corrected(seq_1, alphabet)
            group = 0
            indice = 0

            while group <= len(cluster) - 1 and indice <= len(cluster[group]) - 1:
                name_2 = cluster[group][indice][0]
                pourcentage_id = data_pid[name_1][name_2]
                if pourcentage_id < pid_cluster:
                    group += 1
                    indice = 0
                else:
                    if indice == len(cluster[group]) - 1:
                        cluster[group].append((name_1, seq_1, len_seq_real_1))
                        indice += 2 # avoid infinite loop
                    else:
                        indice += 1
            if group == len(cluster):
                cluster[group] = [(name_1, seq_1, len_seq_real_1)]
    else:
        print(file_seq_non_redondant)
    return cluster



def cluster_representative(cluster):
    """
    Select the first sequence with the longest length in the cluster as the cluster representative
    """
    seq_non_redundant = []
    for group in cluster:
        current_group = cluster[group]
        representative = current_group[0]   
        for elem in current_group:
            if elem[2] > representative[2]: # 2 stands for the corrected lenght of a sequence
                representative = elem
        seq_non_redundant.append(representative[0])   # 0 stands for the name of the sequence
    return seq_non_redundant



def multi_non_redundancy_correction(path_folder_fasta: str, 
                                    path_folder_fasta_non_redondant: str,
                                    path_pid:str,
                                    alphabet: list, pid_cluster: float):
    """Correct the redundancy issue by selecting on representative for each group of sequences
       in a seed with a pid > pid_cluster

    Args:
        path_folder_fasta (str): folder of data to cluster
        path_folder_fasta_non_redondant (str): path where the clustered data is stored
        path_pid (str): path of the pid for each couple of sequences
        alphabet (list): list of character included
        pid_cluster (float): pid min to group the sequences
    """
    folder.creat_folder(path_folder_fasta_non_redondant)

    files = [x for x in Path(path_folder_fasta).iterdir()]
    n_files = len(files)

    start = time.time()
    for file_counter in tqdm(range(n_files), desc = "non-redundant", mininterval=60):
        file = files[file_counter]
        seed = fastaReader.read_multi_fasta(file)

        accession_num = folder.get_accession_number(file)
        path_fasta_non_redondant = f"{path_folder_fasta_non_redondant}/{accession_num}.non_redundant"
        data_pid = np.load(f"{path_pid}/{accession_num}.pid.npy", allow_pickle=True).item()

        cluster = clustering_non_redundant(seed, path_fasta_non_redondant, 
                                           alphabet, pid_cluster, data_pid)
        seq_non_redundant = cluster_representative(cluster)

        with open(file, "r") as file:
            with open(path_fasta_non_redondant, "w") as file_corrected:
                flag_write = False
                for line in file:
                    if line[0] == ">":   
                        if line[1:-1].split(" ")[0] in seq_non_redundant:    # keep the name only
                            flag_write = True
                        else:
                            flag_write = False
                    if flag_write == True:
                        file_corrected.write(line)
    end = time.time()
    print(f"NON-REDUNDANCY: time {'{:_}'.format(round(end - start, 4))} s")


def nbreSeed(path_folder_fasta):
    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()
    count = 0
    for path_file_fasta in files_in_path_folder_fasta:
        len_seed = len(fastaReader.read_multi_fasta(path_file_fasta))
        if len_seed <= 1:
            count += 1       
    print("\nTotal non informative seeds:", count)

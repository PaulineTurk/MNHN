
# IMPORTS
import numpy as np
from tqdm.auto import tqdm
import time
import os

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])

import utils.fastaReader as fastaReader
import utils.folder as folder



# FUNCTIONS

def table_3d_count_initialisation(alphabet):
    table_3d_count = {}
    for aa_origine in alphabet:
        table_3d_count[aa_origine] = {}
        for aa_destination in alphabet:
            table_3d_count[aa_origine][aa_destination] = {}
            for aa_context in alphabet:
                table_3d_count[aa_origine][aa_destination][aa_context]  = 0
    return table_3d_count


# def table_3d_count_origine(seed, nb_seq,
#                            dico_pid_couple, pid_inf, pid_sup,
#                            alphabet,
#                            index_range,
#                            triplet_count, nb_ex_train, delay_num):
#     for i in range(nb_seq-1):
#         name_origine, seq_origine  = seed[i]
#         for j in range(i+1, nb_seq):
#             name_destination, seq_destination  = seed[j]

#             if pid_inf <= dico_pid_couple[name_origine][name_destination] < pid_sup:
#                 for index in index_range:
#                     aa_origine = seq_origine[index]
#                     aa_destination = seq_destination[index]

#                     if all(x in alphabet for x in [aa_origine, aa_destination]):
#                         aa_context = seq_origine[index + delay_num]
#                         if aa_context in alphabet:
#                             nb_ex_train += 1
#                             triplet_count[aa_origine][aa_destination][aa_context] += 1
#                         # cas symétrique
#                         aa_context_sym = seq_destination[index + delay_num]
#                         if aa_context_sym in alphabet:
#                             nb_ex_train += 1
#                             triplet_count[aa_destination][aa_origine][aa_context_sym] += 1

#     return triplet_count, nb_ex_train



def table_3d_count_origine(seed, nb_seq,
                           dico_pid_couple, pid_inf, pid_sup,
                           alphabet,
                           index_range,
                           start, end, step,
                           triplet_count, nb_ex_train):
    for i in range(nb_seq-1):
        name_origine, seq_origine  = seed[i]
        for j in range(i+1, nb_seq):
            name_destination, seq_destination  = seed[j]

            if pid_inf <= dico_pid_couple[name_origine][name_destination] < pid_sup:
                for index in index_range:
                    aa_origine = seq_origine[index]
                    aa_destination = seq_destination[index]

                    if all(x in alphabet for x in [aa_origine, aa_destination]):
                        list_aa_context = [seq_origine[position] for position in range(index + start,
                                                                                   index + end,
                                                                                   step)]
                        if all(x in alphabet for x in list_aa_context):
                            nb_ex_train += 1
                            aa_context = list_aa_context[-1]
                            triplet_count[aa_origine][aa_destination][aa_context] += 1

                        # SYMETRIC CASE
                        list_aa_context_sym = [seq_destination[position] for position in range(index + start,
                                                                                   index + end,
                                                                                   step)]
                        if all(x in alphabet for x in list_aa_context_sym):
                            nb_ex_train += 1
                            aa_context_sym = list_aa_context_sym[-1]
                            triplet_count[aa_destination][aa_origine][aa_context_sym] += 1

    return triplet_count, nb_ex_train



# def table_3d_count_destination(seed, nb_seq,
#                            dico_pid_couple, pid_inf, pid_sup,
#                            alphabet,
#                            index_range,
#                            triplet_count, nb_ex_train, delay_num):
#     """
#     (aa_origine, aa_destination, aa_context) counted with
#     aa_origine, aa_destination and aa_context in alphabet.
#     In addition, aa_origine and aa_destination must be
#     in a pairwise alignment with pid_inf <= pid < pid_sup.
#     destination_right: the contextual amino-acid is at the right
#     of the aa_destination
#     """
#     for i in range(nb_seq-1):
#         name_origine, seq_origine  = seed[i]
#         for j in range(i+1, nb_seq):
#             name_destination, seq_destination  = seed[j]

#             if pid_inf <= dico_pid_couple[name_origine][name_destination] < pid_sup:
#                 for index in index_range:
#                     aa_origine = seq_origine[index]
#                     aa_destination = seq_destination[index]
                      
#                     if all(x in alphabet for x in [aa_origine, aa_destination]):
#                         aa_context = seq_destination[index + delay_num]
#                         if aa_context in alphabet:
#                             nb_ex_train += 1
#                             triplet_count[aa_origine][aa_destination][aa_context] += 1
#                         # cas symétrique
#                         aa_context_sym = seq_origine[index + delay_num]
#                         if aa_context_sym in alphabet:
#                             nb_ex_train += 1
#                             triplet_count[aa_destination][aa_origine][aa_context_sym] += 1

#     return triplet_count, nb_ex_train


def table_3d_count_destination(seed, nb_seq,
                           dico_pid_couple, pid_inf, pid_sup,
                           alphabet,
                           index_range,
                           start, end, step,
                           triplet_count, nb_ex_train):
    for i in range(nb_seq-1):
        name_origine, seq_origine  = seed[i]
        for j in range(i+1, nb_seq):
            name_destination, seq_destination  = seed[j]

            if pid_inf <= dico_pid_couple[name_origine][name_destination] < pid_sup:
                for index in index_range:
                    aa_origine = seq_origine[index]
                    aa_destination = seq_destination[index]

                    if all(x in alphabet for x in [aa_origine, aa_destination]):
                        list_aa_context = [seq_destination[position] for position in range(index + start,
                                                                                   index + end,
                                                                                   step)]
                        if all(x in alphabet for x in list_aa_context):
                            nb_ex_train += 1
                            aa_context = list_aa_context[-1]
                            triplet_count[aa_origine][aa_destination][aa_context] += 1

                        # SYMETRIC CASE
                        list_aa_context_sym = [seq_origine[position] for position in range(index + start,
                                                                                   index + end,
                                                                                   step)]
                        if all(x in alphabet for x in list_aa_context_sym):
                            nb_ex_train += 1
                            aa_context_sym = list_aa_context_sym[-1]
                            triplet_count[aa_destination][aa_origine][aa_context_sym] += 1

    return triplet_count, nb_ex_train



def multi_table_3d_count(path_folder_fasta, path_folder_pid,
                         delay_num, RefSeq, 
                         alphabet, pid_inf, pid_sup,
                         path_NeighborRes):
    """
    Compte le nombre de chaque triplet rencontré dans les alignements multiples (MSA) de path_folder_fasta.
    Seuls les paires d'alignements vérifiant les conditions sur pid et alphabet sont comptés.
    Une unique position de contextage est selectionnée par table de comptage.
    Cette position est déterminée par les paramètres delay_num et kp_SeqChoice.

    La table3d de comptage qui en résulte est enregistrée dans path_NeighborRes
    et est en sortie de cette fonction (ainsi que le nom du dossier des MSA)

    RefSeq: 'o' for origine OR 'd' for destination
    """

    nb_ex_train = 0

    triplet_count = table_3d_count_initialisation(alphabet)

    dico_function = {"o": table_3d_count_origine,
                     "d": table_3d_count_destination}


    files = [x for x in Path(path_folder_fasta).iterdir()]
    nb_files = len(files)

    start = time.time()
    for file_counter in tqdm(range(nb_files), desc='3D count', mininterval=60):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        dico_pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()
        seed = fastaReader.read_multi_fasta(file)
        nb_seq = len(seed)
        if nb_seq > 1:
            len_seq = len(seed[0][1])
            if delay_num < 0:
                index_range = np.arange(- delay_num, len_seq ,1)
                start, end, step = -1, delay_num-1, -1
            else:
                index_range = np.arange(0, len_seq - delay_num, 1)
                start, end, step = 1, delay_num +1, 1

            triplet_count, nb_ex_train = dico_function[RefSeq](seed, nb_seq,
                                                               dico_pid_couple, pid_inf, pid_sup,
                                                               alphabet,
                                                               index_range,
                                                               start, end, step,
                                                               triplet_count, nb_ex_train)
    end = time.time()
    diff = end - start
    items_per_second = nb_files/diff
    print(f"Compte des triplets : {diff:.2f}s | {items_per_second:.2f}it/s")


    # train evaluation
    count_triplet_not_evaluated = 0
    for aa_1 in alphabet:
        for aa_2 in alphabet:
            for aa_3 in alphabet:
                if triplet_count[aa_1][aa_2][aa_3] == 0:
                    count_triplet_not_evaluated += 1
    percentage_triplet_not_evaluated = 100*count_triplet_not_evaluated/(len(alphabet)**3)
    print(f"PARAMETERS NOT ESTIMATED: {percentage_triplet_not_evaluated} %")
    print("N_TRAIN_EXAMPLES :",  '{:_}'.format(nb_ex_train))


    path_triplet_count = f"{path_NeighborRes}/{str(delay_num)}_{RefSeq}"
    np.save(path_triplet_count , triplet_count)

    return triplet_count

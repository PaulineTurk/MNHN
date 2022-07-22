"""
Selection of test examples
"""

# IMPORTS
import numpy as np
import random
import math
from tqdm import tqdm
import time

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.folder as folder


# FUNCTIONS
def example_number_per_seed(path_dico_seed_normal, n_example_test):
    """Pseudo-randomly determine the number of valid examples to pick from each seed

    Args:
        path_dico_seed_normal (str): pre-computed dico storing the proportion of examples
                                     to pick from each seed (to be loaded once)
        n_example_test (int): order of magnitude of the number of valid pairs to pick in
                              total from Pfam_test

    Returns:
        dico: dico storing the number of test examples to pick from each seed
    """
    start = time.time()

    dico_seed_normal = np.load(path_dico_seed_normal, allow_pickle='TRUE').item()
    len_dico_seed_normal = len(dico_seed_normal)

    dico_exemple = {}
    for key in dico_seed_normal:
        dico_exemple[key] = {}
        nb_ex_estimation = dico_seed_normal[key]["WEIGHT_SEED"] * n_example_test
        decimal_part, integer_part = math.modf(nb_ex_estimation)
        proba = random.uniform(0, 1)
        nb_ex_exact = int(integer_part + int(proba < decimal_part))
        dico_exemple[key]["N_Example"] = nb_ex_exact
        dico_exemple[key]["LEN_ALIGN"] = dico_seed_normal[key]["LEN_ALIGN"]

    end = time.time()
    diff = end - start
    items_per_second = len_dico_seed_normal/diff
    print(f'EXAMPLES NUMBER/SEED: time: {diff:.2f}s | {items_per_second:.2f}it/s')

    return dico_exemple


def example_shape(seq_1, seq_2, len_align, position_selected,
                  origine_destination, left_right, context):
    """Formatting of the selected example with completion by *
    in case of overflow of the indices of the alignment

    Args:
        seq_1 (str): first sequence in the pairwise alignment
        seq_2 (str): second sequence in the pairwise alignment
        len_align (float): pairwise alignment number of characters
        position_selected (int): position selected in the pairwise alignment
        origine_destination (str): pick the reference as origine or destination
        left_right (str): pick left or right as the direction of interest
        context (int): max position distance of the contextual amino-acid
                       with the reference

    Returns:
        _type_: _description_
    """
    voisinage = []
    index = position_selected

    if context != 0:
        if origine_destination == "origine":
            seq = seq_1
        if origine_destination == "destination":
            seq = seq_2

        if left_right == "left":
            index -= 1
            while index >= position_selected - context and index in range(len_align):
                voisinage.append(seq[index])
                index -= 1

        else:
            index += 1
            while index <= position_selected + context and index in range(len_align):
                voisinage.append(seq[index])
                index += 1

        out_of_range = context - len(voisinage)
        if out_of_range != 0:
            for i in range(out_of_range):
                voisinage.append("*")

    voisinage = "".join(voisinage)

    return voisinage
       

def random_example_selection(list_example, dico_seq,
                             context_ol, context_or, context_dl, context_dr,
                             nb_ex_test,
                             alphabet):
    """Select of test examples in a seed

    Args:
        list_example (list): list of examples selected
        dico_seq (dico): dictionary of pairwise alignments in a seed
        context_ol (int): index position of the contextual amino-acid
                          with respect to origine/left
        context_or (int): index position of the contextual amino-acid
                          with respect to origine/right
        context_dl (int): index position of the contextual amino-acid
                          with respect to destination/left
        context_dr (int): index position of the contextual amino-acid
                          with respect to destination/right
        nb_ex_test (int): number of examples test to select
        alphabet (list): list of characters of interest

    Returns:
        list: list_example
    """
    # retrieve the list of the pairwise sequence names in a seed
    list_pair_seq_name = tuple(dico_seq.keys())
    # retrive their weight according to "N_VALID_AA_COUPLE"
    weights_list = []
    for key in dico_seq:
        weights_list.append(dico_seq[key]["N_VALID_AA_COUPLE"])
    # select nb_ex_test couples according to their weight
    list_pair_name = random.choices(list_pair_seq_name, weights = weights_list, k = nb_ex_test)

    # select the test examples
    for pair_name in list_pair_name:
        # retrieve the valid indexes
        list_valid_position = dico_seq[pair_name]["LIST_VALID_INDEX"]
        # retrive the pairwise sequences
        if list_valid_position:
            pair_seq = dico_seq[pair_name]["SEQA_SEQB"]
            # random selection of origine/destination for SEQA and SEQB
            proba = random.uniform(0, 1)
            if proba <= 0.5:
                seq_1, seq_2 = pair_seq
            else:
                seq_2, seq_1 = pair_seq
            # random selection of a valid position
            position_selected = random.choice(list_valid_position)
            len_align = len(seq_1)

            # construction of the example and completing the missing positions with *
            voisinage_ol = example_shape(seq_1, seq_2, len_align, position_selected,
                                         "origine", "left", context_ol)
            voisinage_or = example_shape(seq_1, seq_2, len_align, position_selected,
                                         "origine", "right", context_or)
            voisinage_dl = example_shape(seq_1, seq_2, len_align, position_selected,
                                         "destination", "left", context_dl)
            voisinage_dr = example_shape(seq_1, seq_2, len_align, position_selected,
                                         "destination", "right", context_dr)
            # example construction
            example_selected = [seq_1[position_selected],  # aa_1
                                seq_2[position_selected],  # aa_2
                                voisinage_ol,
                                voisinage_or,
                                voisinage_dl,
                                voisinage_dr]
            
            # test de validation de l'exemple:
            example_selected_one_str = "".join(example_selected)
            if all(aa in alphabet for aa in example_selected_one_str):
                list_example.append(example_selected)

    return list_example


def multi_random_example_selection(path_folder_seed, dico_exemple,
                                   path_dico_seq,
                                   context_ol, context_or, context_dl, context_dr,
                                   alphabet):
    """Selection of all the examples for one test

    Args:
        path_folder_seed (str): path of the test seed folder
        dico_exemple (dico): dico with the number of example to pick from each seed
        path_dico_seq (str): path of the folder of pairwise alignments for each seed
        context_ol (int): origine left max index (positive value)
        context_or (int): origine right max index (positive value)
        context_dl (int): destination left max index (positive value)
        context_dr (int): destination right max index (positive value)
        alphabet (list): list of characters of interest

    Returns:
        list: list of the examples for one test
    """
    # initialisation of the list of examples
    list_example = []
    # list of the PosixPath test alignments
    files = [x for x in Path(path_folder_seed).iterdir()]
    nb_files = len(files)
    start = time.time()
    for file_counter in tqdm(range(nb_files), desc='Example selection',
                                   ncols= 100, mininterval=60):
        file = files[file_counter]
        accession_num =  folder.get_accession_number(file)
        if accession_num in dico_exemple:
            nb_ex_test = dico_exemple[accession_num]["N_Example"]
            if nb_ex_test != 0:
                # dico_seq load
                dico_seq = np.load(f"{path_dico_seq}/{accession_num}.seq.npy",
                                   allow_pickle='TRUE').item()
                # selection of nb_ex_test exemples in dico_seq
                list_example = random_example_selection(list_example, dico_seq,
                                                        context_ol, context_or,
                                                        context_dl, context_dr,
                                                        nb_ex_test,
                                                        alphabet)
    end = time.time()
    diff = end-start
    items_per_second = diff/nb_files
    print(f'SELECTION: time: {diff:.2f} s | {items_per_second:.2f} it/s')
    print("N_EXAMPLES:",'{:_}'.format(len(list_example)))

    return list_example

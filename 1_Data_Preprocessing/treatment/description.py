# IMPORTS

import os
import matplotlib.pyplot as plt
import numpy as np
from brokenaxes import brokenaxes

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.fastaReader as fastaReader



#FUNCTIONS

def data_count(path_data, alphabet,
               path_character_percentage,
               path_character_included_percentage):
    """
    path_folder: path of the folder of fasta files to describe
    """
    # initialisation
    n_seed = 0
    n_seq = 0
    n_character = 0
    n_character_included = 0
    character_count = {}

    character_included_count = {}
    for aa in alphabet:
        character_included_count[aa] = 0

    # count
    files = [x for x in Path(path_data).iterdir()]
    nb_files = len(files)
    for file_counter in range(nb_files):
        file = files[file_counter]
        n_seed += 1
        seed = fastaReader.read_multi_fasta(file)
        len_seq = len(seed[0][1])
        for _, seq in seed:
            n_seq += 1
            n_character += len_seq
            for aa in seq:
                if aa in character_count:
                    character_count[aa] += 1
                else:
                    character_count[aa] = 1
                
                if aa in alphabet:
                    n_character_included += 1
                    character_included_count[aa] += 1
                    
    print("\nDESCRIPTION:")
    print(f"N_MULTIPLE ALIGNMENTS: {'{:_}'.format(int(n_seed))}")
    print(f"N_SEQ: {'{:_}'.format(int(n_seq))}")
    print(f"N_CHARACTER: {'{:_}'.format(int(n_character))}")
    print(f"N_CHARACTER_INCLUDED: {'{:_}'.format(int(n_character_included))}")

    print("\nMEANS:")
    # mean len sequence
    if n_seq != 0:
        mean_len_seq = round(n_character_included/n_seq, 2)
        print(f"MEAN_VALID_CHARACTERS_PER_SEQ: {'{:_.2f}'.format(mean_len_seq)}")
    else:
        print("NO SEQUENCE IN THE DATA SELECTED")

    # mean number of sequences per seed
    if n_seed != 0:
        mean_n_seq = round(n_seq/n_seed, 2)
        print(f"MEAN_N_SEQ_PER_SEED: {'{:_.2f}'.format(mean_n_seq)}")
    else:
        print("NO SEED IN THE DATA SELECTED")

    # SAVE
    character_percentage = {k: round(100*v / n_character, 2) for k, v in character_count.items()}
    np.save(path_character_percentage, character_percentage)

    character_included_percentage = {k: round(100*v / n_character_included, 2) for k, v in character_included_count.items()}
    np.save(path_character_included_percentage, character_included_percentage)



def bar_plot_data_description(path_folder_to_describe: str,  path_percentage_description: str, entity_name: str):
    """Description of the distribution of characters in data with bar plot

    Args:
        path_folder_to_describe (str): path of the folder to describe
        path_percentage_description (str): path of the dicitonary
                                           of each character percentage
        entity_name (str): "ACIDE AMINÉ STANDARD" OR "CARACTÈRE"
    """
    data_percentage = np.load(path_percentage_description, allow_pickle=True).item()
    data_percentage_sorted = dict(sorted(data_percentage.items(), key=lambda item: item[1], reverse=True))

    plt.bar(list(data_percentage_sorted.keys()), data_percentage_sorted.values(), color='g')
    plt.xlabel(entity_name)
    plt.ylabel('POURCENTAGE')

    dir_image = os.path.dirname(path_folder_to_describe)
    name_dir = os.path.basename(path_folder_to_describe)

    title_graph = f"Pourcentage de chaque {entity_name} dans \n{name_dir}"
    title_graph_object = f"{dir_image}/{title_graph}"

    plt.title(title_graph)
    plt.savefig(title_graph_object)
    plt.close()

# IMPORTS

import shutil
import time
from sklearn.model_selection import train_test_split

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.folder as folder

# FUNCTIONS

def copy_file_in_folder(file_name, folder_name_source, folder_path_target, extension_file_name_target):
    path_file_source = f"{folder_name_source}/{file_name}"
    accession_num = folder.get_accession_number(path_file_source)

    path_file_target = f"{folder_path_target}/{accession_num}.{extension_file_name_target}"
    shutil.copy2(path_file_source, path_file_target)


def data_split(path_folder_data, path_folder_data_split, percentage_A, name_data_A, name_data_B,
                extension_A = "train", extension_B = "test"):
    """
    Random split files in categorie A and B according to percentage_A

    path_folder_data: folder to split
    path_folder_data_split: folder created where the data splitted are saved
    percentage_A: percentage of files in category A
    name_data_A: name of the folder in path_folder_data_split for the files in category A
    name_data_B: name of the folder in path_folder_data_split for the files in category B
    extension_A: extension of name_data_A
    extension_B: extension of name_data_B
    """
    start = time.time()

    folder.creat_folder(path_folder_data_split)
    path_folder_data_A = f"{path_folder_data_split}/{name_data_A}"
    path_folder_data_B = f"{path_folder_data_split}/{name_data_B}"
    folder.creat_folder(path_folder_data_A)
    folder.creat_folder(path_folder_data_B)

    files = Path(path_folder_data).iterdir()
    data_name = []
    for file_path in files:
        file_name = str(file_path).split("/")[-1]
        data_name.append(file_name)

    fraction_A = percentage_A/100
    files_A, files_B = train_test_split(data_name, train_size = fraction_A)

    for file_name in files_A:
        copy_file_in_folder(file_name, path_folder_data, path_folder_data_A, extension_A)
    
    for file_name in files_B:
        copy_file_in_folder(file_name, path_folder_data, path_folder_data_B, extension_B)

    end = time.time()
    print(f"DATA SPLIT: time {'{:_}'.format(round(end - start, 4))} s")

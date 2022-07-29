# IMPORTS

import sys  
import time
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.folder as folder



# FUNCTION

def multi_capitalization(path_data, path_data_corrected):
    """
    Convert all the lowercase residu into uppercase.

    path_data: path of the folder of fasta file to correct
    path_data_corrected: folder in which the fasta file corrected are saved
    """
    start = time.time()

    folder.creat_folder(path_data_corrected)

    files = [x for x in Path(path_data).iterdir()]
    nb_files = len(files)

    for file_counter in range(nb_files):
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)
        path_file_corrected = f"{path_data_corrected}/{accession_num}.upper"
        with open(file, "r") as file:
            with open(path_file_corrected, "w") as file_corrected:
                for line in file:
                    if line[0] != ">":
                        line = line.upper()
                    file_corrected.write(line)
    
    end = time.time()
    print(f"CAPITALISATION: time {'{:_}'.format(round(end - start, 4))} s")

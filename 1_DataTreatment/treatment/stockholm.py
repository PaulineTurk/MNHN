# IMPORTATIONS

import time
from Bio import AlignIO
from pathlib import Path


import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1])
import utils.folder as folder



# FUNCTIONS


def stockholm_separator(path_file, path_folder_save):
    """
    Separate a multiStockholm file into monoStockholm files.

    path_file_name: path of the multiStockholm file
    path_folder_save: path of th folder where the monoStockholm files generated are saved 
    """
    start = time.time()

    input_handle = open(path_file)
    folder.creat_folder(path_folder_save)

    # collect the accession number of each alignment
    list_accession_num = []
    for l in input_handle:
        if l[0:7] == "#=GF AC":
            init_accession_num = l.index('PF')
            accession_num = l[init_accession_num: -1]
            list_accession_num.append(accession_num)
    input_handle.close()
    nbre_file = len(list_accession_num)

    # generate the monoStockholm files named after the accession number's alignment
    input_handle = open(path_file)
    file_out_nbre = 0
    path_file_out = f"{path_folder_save}/{list_accession_num[file_out_nbre]}.stockholm"
    output_handle = open(path_file_out, "w")

    for l in input_handle:
        output_handle.write(l)
        if l[0:2] == "//" and file_out_nbre <= nbre_file - 2: # avoid generating an empty file at the end
            output_handle.close()
            file_out_nbre += 1
            path_file_out = f"{path_folder_save}/{list_accession_num[file_out_nbre]}.stockholm"
            output_handle = open(path_file_out, "w")

    output_handle.close()
    input_handle.close()

    end = time.time()
    print(f"MULTI TO MONO STOCKHOLM: time {'{:_}'.format(round(end - start, 4))} s")




def multi_stockholm_to_fasta(path_folder_stockholm, path_folder_fasta):
    """
    Convert Stockholm files into Fasta files
    """
    start = time.time()

    folder.creat_folder(path_folder_fasta)

    files_stockholm = Path(path_folder_stockholm).iterdir()
    for file_stockholm in files_stockholm:
        accession_num = folder.get_accession_number(file_stockholm)
        path_file_fasta = f"{path_folder_fasta}/{accession_num}.fasta"
        with open(file_stockholm, "r") as file_stockholm:
            with open(path_file_fasta, "w") as file_fasta:
                alignments = AlignIO.parse(file_stockholm, "stockholm")
                for alignment in alignments:
                    AlignIO.write([alignment], file_fasta, "fasta")

    end = time.time()
    print(f"STOCKHOLM TO FASTA: time {'{:_}'.format(round(end - start, 4))} s")
    
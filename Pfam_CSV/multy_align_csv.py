"""
Write info about pid in a csv file for each Multiple Sequence Alignement
"""
import csv
from pathlib import Path
import os
import time
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta


# global variables

ALPHABET_1 = seq.Alphabet(["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"])

ALPHABET_2 = seq.Alphabet(["B","Z", "X", "O", "U", "-"])

PATH_MSA = "fasta"

# functions

def get_msa(fasta_path):
    """Read an MSA in fasta format

    Args:
        fasta_path (FASTA): MSA in fasta format

    Returns:
        dico: dico of the MSA with sequence name as key and sequence as value
    """
    fasta_file = fasta.FastaFile.read(fasta_path)
    return fasta_file

def get_header(fasta_file):
    """Get the list of the sequences names in an MSA

    Args:
        fasta_file (dico): dictionary of an MSA from get_msa()

    Returns:
        list: name of each sequence in an MSA
    """
    header = dict(fasta_file).keys()
    return header

def pairwise_align(fasta_file, header1, header2):
    """Get a couple of sequences from an MSA

    Args:
        fasta_file (dico): dictionary of an MSA from get_msa()
        header1 (str): a sequence name in the MSA
        header2 (str): a sequence name in the MSA

    Returns:
        tuple: the two string sequences named header1 and header2
               respectively in the MSA
    """
    return fasta_file[header1], fasta_file[header2] # type = str

def pid(sequence1, sequence2, alphabet):
    """Compute the percentage of identity between sequence1 and sequence2
       according to the alphabet chosen

    Args:
        sequence1 (str): one sequence in a MSA
        sequence2 (str): one sequence in a MSA
        alphabet (list): list of the valid characters

    Returns:
        float: percentage of identity
    """
    count_identity = 0
    len_seq_1 = 0
    len_seq_2 = 0
    for char1, char2 in zip(sequence1, sequence2):
        state = 0
        if char1 in alphabet:
            len_seq_1 += 1
            state += 1
        if char2 in alphabet:
            len_seq_2 += 1
            state += 1
        if state == 2 and char1 == char2: # char1 et char2 sont dans ALPHABET_1
            count_identity += 1
    return round(100*(count_identity/min(len_seq_1, len_seq_2)),2)

def get_access_num(path):
    """Get the accession number of a MSA from its path

    Args:
        path (str): path of a MSA

    Returns:
        str: accession number of the MSA
    """
    access_num_1 = os.path.basename(path).split(".")[0]
    access_num_2 = os.path.basename(path).split(".")[1]
    access_num = f"{access_num_1}.{access_num_2}"
    return access_num



# job

files_msa = [x for x in Path(PATH_MSA).iterdir()]
nb_files = len(files_msa)

header_csv = ['access_num', 'name1', 'name2', 'seq1', 'seq2', 'pid']

with open("msa.csv", 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header_csv)

    start = time.time()
    for index, file_msa in enumerate(files_msa):
        accession_number = get_access_num(file_msa)    # nom du MSA
        MSAlign = get_msa(str(file_msa))        # str en arg et non PosixPath
        seq_names = get_header(MSAlign)    # list des noms de seq dans le MSA


        for name1 in seq_names:
            for name2 in seq_names:
                if name1 != name2:
                    seq1, seq2 = pairwise_align(MSAlign, name1, name2)
                    data = [accession_number, name1, name2, seq1, seq2,
                            pid(seq1, seq2, ALPHABET_1)]
                    # write the data
                    writer.writerow(data)
        print(f"{round(100*(index+1)/nb_files, 4)} %")

    end = time.time()
    diff = end - start
    items_per_sec = diff/nb_files
    print(f"csv MSA :{diff:.2f} s | :{items_per_sec:.2f} it/s ")

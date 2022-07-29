import numpy as np
import csv


ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

path_dict = "/home/pauline/Bureau/MNHN_RESULT/3_TABLE_3D/COUNT/6_40_50/or_3.npy"
path_new_csv = "/home/pauline/Bureau/MNHN_RESULT/3_TABLE_3D/COUNT/6_40_50/or_3.csv"

# dict = np.load(path_dict, allow_pickle='True').item()
dict = np.load(path_dict, allow_pickle='True').item()

with open(path_new_csv,'w',
            encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header = ["aa_origin", "aa_destination", "aa_context", "count"]
    writer.writerow(header)

    for aa_o in ALPHABET:
        for aa_d in ALPHABET:
            for aa_c in ALPHABET:
                data = [aa_o, aa_d, aa_c, dict[aa_o][aa_d][aa_c]]
                writer.writerow(data)

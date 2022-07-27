import numpy as np
import csv


ALPHABET = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

# path_dict = "/home/pauline/Bureau/MNHN_RESULT/3_TABLE_3D/COUNT_40_50_6/ol_3.npy"
path_dict_COMPARE = "/home/pauline/Bureau/MNHN_RESULT/3_TABLE_3D/COUNT_40_50_6_test_concat/or_1.npy"
path_new_csv = "/home/pauline/Bureau/MNHN_RESULT/3_TABLE_3D/COUNT_40_50_6_test_concat/or_1_new.csv"

# dict = np.load(path_dict, allow_pickle='True').item()
dict_COMPARE = np.load(path_dict_COMPARE, allow_pickle='True').item()

with open(path_new_csv,'w',
            encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header = ["aa_origin", "aa_destination", "aa_context"]
    writer.writerow(header)

    for aa_o in ALPHABET:
        for aa_d in ALPHABET:
            for aa_c in ALPHABET:
                data = [aa_o, aa_d, aa_c, dict_COMPARE[aa_o][aa_d][aa_c]]  #  - dict_COMPARE[aa_o][aa_d][aa_c]
                writer.writerow(data)

"""
Selection des exemples test/train et conversion des lettre en numériques
"""

import csv
import pandas as pd

CODE = {"A":1, "R":2, "N":3, "D":4, "C":5,
        "Q":6, "E":7, "G":8, "H":9, "I":10,
        "L":11, "K":12, "M":13, "F":14, "P":15,
        "S":16, "T":17, "W":18, "Y":19, "V":20,
        "B":21,"Z":22, "X":23, "O":24, "U":25,
        "-": 26,
        "*": 27}

def encoding(seq, code):
    """Encode a sequence according to the code transformation

    Args:
        seq (str): a sequence in a MSA
        code (dico): dictionary to code each letter in the sequence
                    into a numerical value

    Returns:
        list: a list of floats i.e the numeric equivalent of the initial
              sequence
    """
    list_num = []
    for char in seq:
        list_num.append(int(code[char]))
    return list_num




################################################################################


PATH_MSA = "msa.csv"

header_csv = ['access_num', 'name1', 'name2', 'seq1', 'seq2', 'pid']

col_names = pd.read_csv(PATH_MSA, nrows=0).columns
types_dict = {'access_num': str, 'name1': str, 'name2': str,
              'seq1': str, 'seq2':str, 'pid':float}
types_dict.update({col: str for col in col_names if col not in types_dict})
data = pd.read_csv(PATH_MSA, delimiter = ",", dtype=types_dict)

# print(data)


PID_INF = 50  # ajouter avec arg.parse
PID_SUP = 100
NUM_MAX = 20

sub_data = data.loc[(data['pid'] >= PID_INF) & (data['pid'] < PID_SUP)]
# print(sub_data)

header_csv = ['num1', 'num2', 'pid']

with open(f"context_free_{PID_INF}_{PID_SUP}.csv",
        'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header_csv)

    for i in sub_data.index:
        seq1 = encoding(sub_data["seq1"][i], CODE)
        seq2 = encoding(sub_data["seq2"][i], CODE)

        for num1, num2 in zip(seq1, seq2):
            if int(num1) <= 20 and int(num2) <= 20:
                data = [num1, num2, sub_data["pid"][i]]
                writer.writerow(data)

# nb_lines in csv
with open(f"context_free_{PID_INF}_{PID_SUP}.csv", encoding='UTF8') as f:
    reader = csv.reader(f)
    list_reader = list(reader)
    print(f"nb ex : {len(list_reader)}")

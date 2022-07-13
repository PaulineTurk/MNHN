"""
Creat a .csv with destination amino-acid is the max in digit code:
python3 max_neighbour.py
"""

import csv
import random

N_EXAMPLES = 45_000
N_NEIGHBOUR = 3
LIST_CHARACTERS = [0, 1, 2, 3, 4, 5, 6, 7, 8]
NAME_CSV = "max_neighbour"
with open(f"{NAME_CSV}.csv",
    'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    for i in range(N_EXAMPLES):
        data = [random.choice(LIST_CHARACTERS)]
        for j in range(N_NEIGHBOUR):
            data.append(random.choice(LIST_CHARACTERS))
        data.append(max(data[1:]))
        data = tuple(data)
        writer.writerow(data)

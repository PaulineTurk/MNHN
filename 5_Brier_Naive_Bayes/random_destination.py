"""
Creat a .csv with uniform random destination amino-acid:
python3 random_destination.py
"""

import csv
import random

N_EXAMPLES = 45000
LIST_CHARACTERS = [0, 1, 2, 3, 4, 5, 6, 7, 8]
N_EXAMPLES_PER_CHARACTER = int(N_EXAMPLES/len(LIST_CHARACTERS))

NAME_CSV = "random_destination"
with open(f"{NAME_CSV}.csv",
    'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    for num in LIST_CHARACTERS:
        for i in range(N_EXAMPLES_PER_CHARACTER):
            data = num, random.choice(LIST_CHARACTERS)
            writer.writerow(data)

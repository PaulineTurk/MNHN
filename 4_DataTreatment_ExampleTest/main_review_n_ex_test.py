"""
REVIEW: Preprocessing of data test selection
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

# latex activation
rc('text', usetex=True)

import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])

# PARAMETERS
DATA = f"{file.parents[2]}/MNHN_RESULT/4_DATA_EXEMPLE_TEST/REVIEW"
NAME_CSV_FILE = "ex_test.csv"

print("_______________________________________________________________________")
print("                          REVIEW PRE-PROCESSING                        ")
print("_______________________________________________________________________")
# data loading
data_csv = np.loadtxt(f"{DATA}/{NAME_CSV_FILE}" , delimiter=",", dtype=np.float32, skiprows=0)

# conversion to pandas dataframe + headers addition
df = pd.read_csv(f"{DATA}/{NAME_CSV_FILE}",
                 names=['pid_inf', 'pid_sup', 'n_valid_seed', 'n_non_info_seed', 'n_valid_aa_couple_global'])

# sort by increasing pid_inf
df = df.sort_values('pid_inf')
print (df)

# y-axis in bold
rc('font', weight='bold')

# Values of each group in percentage
n_aa_couple_no_context = sum([elem for elem in df['n_valid_aa_couple_global']])
print(n_aa_couple_no_context)
bars_n_aa_couple = [100*elem/n_aa_couple_no_context for elem in df['n_valid_aa_couple_global']]
print(bars_n_aa_couple)

 
# The position of the bars on the x-axis
n_lines = len(bars_n_aa_couple)
r = [3*i for i in range(n_lines)]
 
# Names of group and bar width
names = []
for pid_inf, pid_sup in zip(df['pid_inf'], df['pid_sup']):
    names.append(f"{int(pid_inf)}-{int(pid_sup)}")
print(names)
barWidth = 3


plt.figure(figsize=(8,5))

# Create bars
graph = plt.bar(r, bars_n_aa_couple, edgecolor='white',
                width=barWidth, alpha=0.7)

# Custom axis + legend + title
plt.xticks(r, names, fontweight='bold')
plt.xlabel("Intervalle de PID", fontsize=13)
plt.ylabel(f"\% d'exemples tests non-contextuels", fontsize=13)
n_aa_couple_no_context = '{:,}'.format(n_aa_couple_no_context).replace(',', ' ')
title = f"Répartition du pourcentage d'exemples tests non-contextuels\n selon l'intervalle de PID (TOTAL: {n_aa_couple_no_context} exemples tests non-contextuels)"
plt.title(title, loc='center', fontsize=16)

i = 0
for p in graph:
    width = p.get_width()
    height = p.get_height()
    x, y = p.get_xy()
    plt.text(x+width/2,
             y+height*1.01,
             str(round(bars_n_aa_couple[i],2))+'\%',
             ha='center',
             weight='bold')
    i+=1

# plt.show()
# Save graphic
plt.savefig(f"{DATA}/graph_percentage_msa.png")

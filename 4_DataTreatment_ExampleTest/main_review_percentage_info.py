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
n_seed = df['n_valid_seed'][0] + df['n_non_info_seed'][0]
print(n_seed)
bars_non_info_seed = [100*elem/n_seed for elem in df['n_non_info_seed']]
print(bars_non_info_seed)
bars_valid_seed = [100*elem/n_seed for elem in df['n_valid_seed']]
print(bars_valid_seed)
 
# The position of the bars on the x-axis
n_lines = len(bars_non_info_seed)
r = [3*i for i in range(n_lines)]
 
# Names of group and bar width
names = []
for pid_inf, pid_sup in zip(df['pid_inf'], df['pid_sup']):
    names.append(f"{int(pid_inf)}-{int(pid_sup)}")
print(names)
barWidth = 3


plt.figure(figsize=(8,5))

# Create bars (inf)
plt.bar(r, bars_non_info_seed, edgecolor='white', width=barWidth, hatch = '/',
        label = 'NON-INFORMATIF', color = "tab:red", alpha=0.7) #, color='#7f6d5f'
# Create bars (sup)
plt.bar(r, bars_valid_seed, bottom=bars_non_info_seed, edgecolor='white', width=barWidth,
        label = 'INFORMATIF', color = "tab:green", alpha=0.7) # ,color='#557f2d'

# Custom axis + legend + title
plt.xticks(r, names, fontweight='bold')
plt.xlabel("Intervalle de PID", fontsize=13)
plt.ylabel(f"\% d'alignements multiples", fontsize=13)
n_seed = '{:,}'.format(n_seed).replace(',', ' ')
title = f"Répartition de l'origine des exemples tests\n selon l'intervalle de PID (TOTAL: {n_seed} alignements multiples)"
plt.title(title, loc='center', fontsize=18)
plt.legend()  # loc=9

# plt.show()
# Save graphic
plt.savefig(f"{DATA}/graph_percentage_info_msa.png")

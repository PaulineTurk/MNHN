"""
python3 hist_pc_3d.py
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


file = Path(__file__).resolve()
sys.path.append(file.parents[0])




LIST_PSEUDO_COUNTER_3D = [0,
                          pow(10, -3),
                          pow(10, -2),
                          pow(10, -1),
                          pow(10, 0),
                          pow(10, 1),
                          pow(10, 2),
                          pow(10, 3)]


EXP_NAME = "5_PC_3D_SELECTION_CORRECTION_BIAIS_test_3"

MAIN_DATA = f"{file.parents[2]}/MNHN_RESULT/{EXP_NAME}/Experiment_Brier_PC_3D_SELECTION_60_70"
# FOLDER REGISTRATION
path_image = f"{MAIN_DATA}/REVIEW"
if not os.path.exists(path_image):
    os.makedirs(path_image)



f,a = plt.subplots(3,6)
a = a.ravel()
index = 0
for ax in a:
    if index < len(LIST_PSEUDO_COUNTER_3D):
        PC = float(LIST_PSEUDO_COUNTER_3D[index])
        DATA = f"{MAIN_DATA}/convergence_{PC}.csv"

        df = pd.read_csv(DATA,
                    names=['context', 'exemple', 'brier_unit', 'brier_unit_no_c'],
                    dtype=str)
        
        x = [float(elem) for elem in df['brier_unit'][1:].tolist()]
        y = [float(elem) for elem in df['brier_unit_no_c'][1:].tolist()]
        bins = np.linspace(0, 2, 100)

        ax.hist(x, bins, alpha=0.5, label='C')
        ax.hist(y, bins, alpha=0.5, label='/C')
        ax.legend(loc='upper right')
        ax.set_title(f"{PC}")

        index += 1



f.delaxes(ax)


plt.show()

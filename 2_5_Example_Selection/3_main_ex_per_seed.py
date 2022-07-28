import shutil
import argparse

import pandas as pd

# PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("path_original",
    help="path of the initial csv of examples per seed", type=str)
parser.add_argument("path_target",
    help="path of the new csv of fraction of examples per seed", type=str)
args = parser.parse_args()

path_original = r'/home/pauline/Bureau/MNHN_RESULT/2_5_EXAMPLES/num_ex.csv'
path_target = r'/home/pauline/Bureau/MNHN_RESULT/2_5_EXAMPLES/frac_ex.csv'

shutil.copyfile(path_original, path_target)

df = pd.read_csv (path_target)

total_num_ex = sum(df["num_ex"])
print(f"N_EX_TOTAL: {total_num_ex}")
df["frac_ex"] = df["num_ex"]/total_num_ex
df.to_csv(path_target, index=False)

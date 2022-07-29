# IMPORTS

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


shutil.copyfile(args.path_original, args.path_target)

df = pd.read_csv (args.path_target)

total_num_ex = sum(df["num_ex"])
print(f"N_EX_TOTAL: {total_num_ex}")
df["frac_ex"] = df["num_ex"]/total_num_ex
df.to_csv(args.path_target, index=False)

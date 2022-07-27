import csv
import shutil

import pandas as pd


path_original = r'/home/pauline/Bureau/MNHN_RESULT/2_5_EXAMPLES/num_ex.csv'
path_target = r'/home/pauline/Bureau/MNHN_RESULT/2_5_EXAMPLES/frac_ex.csv'

shutil.copyfile(path_original, path_target)

df = pd.read_csv (path_target)

# remove the lines with no valid examples
df = df.drop(df[df.num_ex == 0].index)

total_num_ex = sum(df["num_ex"])
print(total_num_ex)  # 672_510_626
df["frac_ex"] = df["num_ex"]/total_num_ex
df.to_csv(path_target, index=False)

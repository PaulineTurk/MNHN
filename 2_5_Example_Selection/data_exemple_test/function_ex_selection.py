import shutil
import os

import pandas as pd



def ex_selection_ol(original, pid_inf, pid_sup,
                    target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    df_restricted = df.loc[(df['pid'] >= pid_inf) & (df['pid'] < pid_sup)]
    columns_ol = ["pid", "name_origin", "name_destination"] + [f"aa_ol_{i}" for i in range(1, L+1)]
    df_restricted_ol = df_restricted[columns_ol]
    df_restricted_ol.to_csv(target, index=False)
    
    
def ex_selection_or(original, pid_inf, pid_sup,
                    target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    df_restricted = df.loc[(df['pid'] >= pid_inf) & (df['pid'] < pid_sup)]
    columns_or = ["pid", "name_origin", "name_destination"] + [f"aa_or_{i}" for i in range(1, L+1)]
    df_restricted_or = df_restricted[columns_or]
    df_restricted_or.to_csv(target, index=False)

def ex_selection_dl(original, pid_inf, pid_sup,
                    target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    df_restricted = df.loc[(df['pid'] >= pid_inf) & (df['pid'] < pid_sup)]
    columns_dl = ["pid", "name_origin", "name_destination"] + [f"aa_dl_{i}" for i in range(1, L+1)]
    df_restricted_dl = df_restricted[columns_dl]
    df_restricted_dl.to_csv(target, index=False)

def ex_selection_dr(original, pid_inf, pid_sup,
                    target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    df_restricted = df.loc[(df['pid'] >= pid_inf) & (df['pid'] < pid_sup)]
    columns_dr = ["pid", "name_origin", "name_destination"] + [f"aa_dr_{i}" for i in range(1, L+1)]
    df_restricted_dr = df_restricted[columns_dr]
    df_restricted_dr.to_csv(target, index=False)


# original = r'/Users/pauline/Desktop/MNHN_RESULT_MINI/2_5_EXAMPLES/EXAMPLES_6/PF09847.12.csv'
# path_target = "/Users/pauline/Desktop/MNHN_RESULT_MINI/2_5_EXAMPLES/EXAMPLES_6_restricted_dr"
# os.makedirs(path_target, exist_ok=True)

# target = f'{path_target}/PF09847.12.csv'
# pid_inf, pid_sup = 30, 40
# L = 6


# # ex_selection_or(original, pid_inf, pid_sup, target, L)
# # ex_selection_ol(original, pid_inf, pid_sup, target, L)
# # ex_selection_dl(original, pid_inf, pid_sup, target, L)
# ex_selection_dr(original, pid_inf, pid_sup, target, L)
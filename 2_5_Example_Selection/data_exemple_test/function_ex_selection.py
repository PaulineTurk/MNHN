import shutil
import os

import pandas as pd

# ouverture fichier in main

def ex_selection_ol(original, target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    columns_ol = ["aa_origin", "aa_destination"] + [f"aa_ol_{i}" for i in range(1, L+1)]
    df_restricted_ol = df[columns_ol]
    df_restricted_ol.to_csv(target, index=False)
    
def ex_selection_or(original, target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    columns_or = ["aa_origin", "aa_destination"] + [f"aa_or_{i}" for i in range(1, L+1)]
    df_restricted_or = df[columns_or]
    df_restricted_or.to_csv(target, index=False)


def ex_selection_dl(original, target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    columns_dl = ["aa_origin", "aa_destination"] + [f"aa_dl_{i}" for i in range(1, L+1)]
    df_restricted_dl = df[columns_dl]
    df_restricted_dl.to_csv(target, index=False)
    
def ex_selection_dr(original, target, L):
    shutil.copyfile(original, target)
    df = pd.read_csv(target)
    columns_dr = ["aa_origin", "aa_destination"] + [f"aa_dr_{i}" for i in range(1, L+1)]
    df_restricted_dr = df[columns_dr]
    df_restricted_dr.to_csv(target, index=False)
    
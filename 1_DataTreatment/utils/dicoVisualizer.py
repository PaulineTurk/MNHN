import numpy as np
import pandas as pd

def dico_visualizer(path_dico):
    """View and import a dictionary from its path

    Args:
        path_dico (str): path of the dico

    Returns:
        dico: the dictionary of interest
    """
    dico = np.load(path_dico, allow_pickle='TRUE').item()
    df_dico = np.transpose(pd.DataFrame.from_dict(dico))
    print(path_dico)
    print(df_dico)
    return df_dico

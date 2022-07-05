import numpy as np
from statsmodels.stats.weightstats import ztest as ztest


list_standard_aa       = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

# test doit-etre fait sur table3d compte? 
# ou table3d proba? ou pas d'importance?

# ex test des tables3d proba m1
path_tables_3d = "/home/pauline/Bureau/MNHN_RESULT/3_TABLE_3D_M1"


def calculArrayDiffkp(path_tables_3d, position, pid_inf, pid_sup, list_standard_aa):
    """
    position: décalage de position, négatif si à gauche, positif si à droite
    """

    table3d_k = np.load(f"{path_tables_3d}/table_3d_{pid_inf}_{pid_sup}/proba_cond_({position},k).npy", allow_pickle='TRUE').item()
    table3d_p = np.load(f"{path_tables_3d}/table_3d_{pid_inf}_{pid_sup}/proba_cond_({position},p).npy", allow_pickle='TRUE').item()

    len_list_standard_aa = len(list_standard_aa)

    # intialisation de l'array de la différence terme à terme table3d_k - table3d_p
    array_diff_k_p = np.empty(len_list_standard_aa**3, dtype=object)

    # calcul de array_diff_k_p
    i = 0
    for aa_1 in list_standard_aa:
        for aa_2 in list_standard_aa:
            for aa_c in list_standard_aa:
                array_diff_k_p[i] = table3d_k[aa_1][aa_2][aa_c] - table3d_p[aa_1][aa_2][aa_c]
                i += 1

    #print("array_diff_k_p: ", array_diff_k_p)
    return array_diff_k_p


def ztestEval(array_diff_k_p):
    mean_diff = np.mean(array_diff_k_p)
    ztest_Score, p_value = ztest(array_diff_k_p, 
                                 value = mean_diff,
                                 alternative='two-sided')
    if p_value != 1:
        print("mean_diff :", mean_diff)
        print("ztest_Score :", ztest_Score)
        print("p_value :", p_value)


for j in range(-10,11):
    if j!=0:  # pas de distinction o/d dans la cas des tables 2D
        for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
            position, pid_inf, pid_sup = j, i, i + 10
            print(j, i, i +10)
            array_diff_k_p = calculArrayDiffkp(path_tables_3d, position, pid_inf, pid_sup, list_standard_aa)
            ztestEval(array_diff_k_p)

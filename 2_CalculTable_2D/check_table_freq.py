import numpy as np

path = "/home/pauline/Bureau/MNHN_RESULT/2_TABLE_2D/90_100_1/count.npy"

data = np.load(path, allow_pickle=True).item()

print(data["A"]["L"])

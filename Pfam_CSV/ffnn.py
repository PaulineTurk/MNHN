"""
FFNN
"""

import torch
import torch.nn.functional as F
import numpy as np

NB_NEIGHBOUR = 3

data_csv = np.loadtxt("context_free_50_100.csv", delimiter=",",
 dtype=np.float32, skiprows=1)
# pylint: disable=E1101
data_tensor = torch.from_numpy(data_csv)
# pylint: enable=E1101

data_tensor = data_tensor.long()
pred = F.one_hot(data_tensor[:, [0,2]], 21)
target = F.one_hot(data_tensor[:, 1], 21)
print(pred)
print(target)

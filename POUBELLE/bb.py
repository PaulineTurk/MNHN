"""
FFNN _ bb
"""

import math
import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import torch.utils.data as data_utils
import numpy as np
from sklearn.model_selection import train_test_split


m = nn.Softmax(dim=1)
input = torch.randn(2, 3)
output = m(input)
print(output)

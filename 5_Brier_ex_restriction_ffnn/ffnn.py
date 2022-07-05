"""
FFNN
"""

import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader # il s'agit de 2 classes
import torch.utils.data as data_utils
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

N_NEIGHBOUR = 3
N_CHARACTER = 20


#######################################


# device config
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# hyper parameters
# INPUT_SIZE = N_NEIGHBOUR + 1
HIDDEN_SIZE = 100   # comment le choisir ?
N_CLASSES = N_CHARACTER
N_EPOCHS = 2
BATCH_SIZE = 3
INPUT_SIZE = (N_CHARACTER+1)*(N_NEIGHBOUR+1)*BATCH_SIZE # ? pas du tout sure de ca ...
LEARNING_RATE = 0.001
TEST_SIZE = 0.2
RANDOM_STATE = 1234 # ne pas le fixer ?
DATA = "data_70_80_multi_origine_gauche.csv"


################################################################
# DATASET
class AminoAcidDataset(Dataset):
    """Construction of the AminoAcide Dataset of examples

    Args:
        Dataset (csv): file with an example per line
    """

    def __init__(self, data_path:str, n_classes:int):
        # data loading
        data_csv = np.loadtxt(data_path, delimiter=",", dtype=np.float32, skiprows=0)
        # skiprows=0, csv with no headers

        # pylint: disable=E1101
        data_tensor = torch.from_numpy(data_csv) # numpy to tensor
        # pylint: enable=E1101
        data_tensor = data_tensor.long() # int68
        self.n_samples = data_tensor.shape[0]
        self.n_columns = data_tensor.shape[1]
        self.n_features = self.n_columns-1 # remove the target column

        # slicing des données selon X et Y
        # avec OneHot embedding
        column_x = [i for i in range(self.n_columns) if i != 1]
        self.X = F.one_hot(data_tensor[:, column_x], n_classes +1) # compte commence à 0
                                                           # peut-etre qu'il faudrait
                                                           # ne pas avoir besoin de ca ...
                                                            
        # à automatiser le nombre de caractères ...?
        self.Y = F.one_hot(data_tensor[:, 1], n_classes +1)

    
    def __getitem__(self, index):
        # dataset[index]
        return self.X[index], self.Y[index]    # chaque ligne c un tuple: x, y   i.e features, label

    def __len__(self):
        # len(dataset)
        return self.n_samples


# création de l'instance du dataset
dataset = AminoAcidDataset(DATA, N_CLASSES)
# print(data.n_samples, data.n_features)

X_train, X_test, Y_train, Y_test = train_test_split(dataset.X,
                                                    dataset.Y,
                                                    test_size=TEST_SIZE,
                                                    random_state=RANDOM_STATE)

first_data = dataset[0] # premier ex récupérer grace à la méthode __getitem__
X_first_data, Y_first_data = first_data
# print(f"nbre de features dans dataset: {dataset.n_features}")



################################################################
# DATALOAD

# création de l'instance d'objet dataloader

################# c faux ..............

# print(X_train.shape)
# print(Y_train.shape)
# train_data = torch.hstack((X_train, Y_train))
# print(train_data.shape)

dataset_train = data_utils.TensorDataset(X_train, Y_train)
dataloader_train = DataLoader(dataset=dataset_train, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)
# num_workers: pour utiliser plusieurs processeurs et donc accélérer le dataloading ...
# jusqu'à combien je peux monter?


dataset_test = data_utils.TensorDataset(X_test, Y_test)
dataloader_test = DataLoader(dataset=dataset_test, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)



# dataloader_train_iter = iter(dataloader_train)  # peut-etre pas besoin de passer en mode .iter() ?
# print(dataloader_train_iter.next())

# utilisation de l'objet dataloader
# dataiter = iter(dataloader)
# data = dataiter.next()  # take the next item
# print(data)
# X_data, Y_data = data
# print(X_data, Y_data)


# train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=batch_size,
#                                             shuffle=True)

# test_loader = torch.utils.data.DataLoader(dataset=test_dataset, batch_size=batch_size,
#                                             shuffle=False)  # inutil de shuffle pour l'évaluation









################################################################
# STRUCTURE SIMPLE AVEC UNE HIDDEN LAYER
class NeuralNet(nn.Module):

    def __init__(self, input_size, hidden_size, n_classes):
        super(NeuralNet, self).__init__()   # hérite d'elle-meme ?

        # initialisation de chaque layer
        self.l1 = nn.Linear(input_size, hidden_size)   # input_size, output_size
        print(input_size, hidden_size)
        self.relu = nn.ReLU() # fonction d'activation, pour la hidden layer (c quoi sa logique?)
        self.l2 = nn.Linear(hidden_size, n_classes)   # input_size, output_size
        print(hidden_size, n_classes)

    def forward(self, x):  # x est l'input
        out = self.l1(float(x.items())) # car long imposé pour la conversion en oneHot
                                # alors qu'un FLoat est attendu dans cette méthode
        out = self.relu(out) # prend en entrée la sortie précédente
        out = self.l2(out) # prend en entrée la sortie précédente

        # pas la classique softmax à la fin, on veut appliquer Brier au vecteur de proba
        return out


## modèle
model = NeuralNet(INPUT_SIZE, HIDDEN_SIZE, N_CLASSES).to(device)  # send the model to the device

### LOSS
# criterion = nn.CrossEntropyLoss() # le softmax est implémenté dedans
criterion = nn.MSELoss()  # Mean Square Error --> prob de regression
# criterion = nn.NLLLoss() # Negative Log Likelihood --> prob de classification
# optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)

### OPTIMIZER: Stochastic Gradient Descent
optimizer = torch.optim.SGD(model.parameters(), lr=LEARNING_RATE)


## training loop
n_total_steps = len(dataloader_train)  # nb de batch par epoch # histoire de math.ceil ...?
for epoch in range(N_EPOCHS):
    for i, (X, Y) in enumerate(dataloader_train):
        # reshape input from 100, 1, 28, 28 to 100, 784
        # X = X.reshape(-1, 28*28).to(device)   # je ne sais pas si le reshape est obligatoire ?
                                                # si oui, à le faire en prenant les bonnes dimensions
        X = X.reshape(-1, 12*21).to(device) # l'envoie à device en corrigeant
                                            # si nécessaire le format
        Y = Y.to(device)

        # forward
        Y_pred = model(X)
        loss = criterion(Y_pred, Y)

        # backward
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (i+1) % 100 == 0:
            print(f"epoch {epoch+1} / {N_EPOCHS}, step {i+1}/{n_total_steps}, loss = {loss.item():.4f}")

    

# # test
# with torch.no_grad():  # testing is very fast
#     n_correct = 0
#     n_samples = 0

#     for images, labels in test_loader:
#         images = images.reshape(-1, 28*28).to(device)
#         labels = labels.to(device)
#         outputs = model(images)  #calcul de la prédiction sur le modèle entrainé


#         # dans mon cas je veux les outputs pour calculer brier dessus

#         # value, index
#         _, predictions = torch.max(outputs, 1)  # predictions: classes des prédictions faite
#                                                 # pour chaque image test
#         n_samples += labels.shape[0]  # nb of axamples in the current batch, should be 100 here
#         n_correct += (predictions == labels).sum().item() # .item(): pour prendre l'unique valeur du tenseur 

#     acc = 100.0 * n_correct / n_samples
#     print(f"accuracy: {acc}")

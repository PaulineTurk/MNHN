"""
FFNN: python3 ffnn.py random_destination.csv > random_destination_$$.txt 2>&1
"""
import argparse
import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import torch.utils.data as data_utils
import numpy as np
from sklearn.model_selection import train_test_split

parser = argparse.ArgumentParser()
parser.add_argument("path_data", help="dataset global path", type=str)
args = parser.parse_args()

N_NEIGHBOUR = 0
N_CHARACTER = 9   # pour le test


#######################################

print("________________________________________________")
print("             LEARNING IDENTITY TEST             ")
print("           aa_origine == aa_destination         ")
print("________________________________________________")

# device config
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')




# hyper parameters
N_CLASSES = N_CHARACTER
INPUT_SIZE = N_CHARACTER # ? pas du tout sure de ca ...
HIDDEN_SIZE = 100   # comment le choisir ?
N_EPOCHS = 100
BATCH_SIZE = 50
LEARNING_RATE = 0.001
TEST_SIZE = 0.2
RANDOM_STATE = 1234 # ne pas le fixer ?
DATA = args.path_data






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
        self.X = F.one_hot(data_tensor[:, column_x], n_classes) # compte commence à 0
                                                           # peut-etre qu'il faudrait
                                                           # ne pas avoir besoin de ca ...                                                     
        # à automatiser le nombre de caractères ...?
        self.Y = F.one_hot(data_tensor[:, 1], n_classes)

    
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

dataset_train = data_utils.TensorDataset(X_train, Y_train)
dataloader_train = DataLoader(dataset=dataset_train, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)
# num_workers: pour utiliser plusieurs processeurs et donc accélérer le dataloading ...
# print(iter(dataloader_train).next())
# print(iter(dataloader_train).next())
# print(iter(dataloader_train).next())


dataset_test = data_utils.TensorDataset(X_test, Y_test)
dataloader_test = DataLoader(dataset=dataset_test, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)
print("")
print("____________________")
print("  DATASET/DATALOAD  ")
print("____________________")
print("")
print("DATASET with 9 amino-acids encoded from 0 to 8")
print(f"N_NEIGHBOUR : {N_NEIGHBOUR}")
print(f"N_EXAMPLES (TOTAL): {len(dataset)}")
print(f"Example of one X_train, Y_train: {first_data}")
print(f"N_EXAMPLES_TRAIN: {len(dataset_train)}")
print(f"N_EXAMPLES_TEST: {len(dataset_test)}")
print("(aa_origine == aa_destination) == True  # for all examples")

print("")
print("____________________")
print("  HYPER-PARAMETERS  ")
print("____________________")
print("")
print(f"N_CLASSES : {N_CLASSES}")
print(f"INPUT_SIZE : {INPUT_SIZE}")
print(f"HIDDEN_SIZE : {HIDDEN_SIZE}")
print(f"N_EPOCHS : {N_EPOCHS}")
print(f"BATCH_SIZE : {BATCH_SIZE}")
print(f"LEARNING_RATE : {LEARNING_RATE}")
print(f"TEST_SIZE : {TEST_SIZE}")
print(f"RANDOM_STATE : {RANDOM_STATE}")

################################################################
# STRUCTURE SIMPLE AVEC UNE HIDDEN LAYER
class NeuralNet(nn.Module):

    def __init__(self, input_size, hidden_size, n_classes):
        super(NeuralNet, self).__init__()   # hérite d'elle-meme ?

        # initialisation de chaque layer
        self.l1 = nn.Linear(input_size, hidden_size)   # input_size, output_size
        # print(input_size, hidden_size)
        self.relu = nn.ReLU() # fonction d'activation, pour la hidden layer (c quoi sa logique?)
        self.l2 = nn.Linear(hidden_size, n_classes)   # input_size, output_size
        # print(hidden_size, n_classes)
        self.m = nn.Softmax(dim=1)

    def forward(self, x):  # x : RNN input
        out = self.l1(x)
        out = self.relu(out)
        out = self.l2(out)
        out = self.m(out)

        # pas la classique softmax à la fin, on veut appliquer Brier au vecteur de proba
        return out


## modèle
model = NeuralNet(INPUT_SIZE, HIDDEN_SIZE, N_CLASSES).to(device)  # send the model to the device

### LOSS
# criterion = nn.CrossEntropyLoss() # le softmax est implémenté dedans
criterion = nn.MSELoss()  # Mean Square Error --> prob de regression
# criterion = nn.NLLLoss() # Negative Log Likelihood --> prob de classification


### OPTIMIZER: Stochastic Gradient Descent
#optimizer = torch.optim.SGD(model.parameters(), lr=LEARNING_RATE)
optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
print("")
print("____________________")
print("        MODEL       ")
print("____________________")
print("")
print("Feed Forward Neural Network with 1 HIDDEN LAYER + a final SOFTMAX LAYER")
print("LOSS: MSELoss")
print("OPTIMIZER: Adam")




## training loop
print("")
print("____________________")
print("   TRAINING LOOP    ")
print("____________________")
print("")
total_samples = len(dataloader_train)  # nb de batch par epoch # histoire de math.ceil ...?
for epoch in range(N_EPOCHS):
    for i, (X, Y) in enumerate(dataloader_train):

        # reshape input from 100, 1, 28, 28 to 100, 784
        # X = X.reshape(-1, 28*28).to(device)   # je ne sais pas si le reshape est obligatoire ?
                                                # si oui, à le faire en prenant les bonnes dimensions
        X = X.reshape(-1, INPUT_SIZE).to(device) # l'envoie à device en corrigeant
                                            # si nécessaire le format
        Y = Y.to(device)
        # print(f"epoch : {epoch}")
        # print(f"X : \n{X}")
        # print(f"Y : \n{Y}")

        # forward
        Y_pred = model(X.float())
        loss = criterion(Y_pred, Y.float())

        # backward
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (i+1) % 100 == 0:
            print(f"epoch {epoch+1} / {N_EPOCHS}, step {i+1}/{total_samples},\
                    loss = {loss.item():.4f}")

    

# # test
print("")
print("____________________")
print("    TESTING LOOP    ")
print("____________________")
print("")
with torch.no_grad():  # testing is very fast
    N_CORRECT = 0
    N_SAMPLES = 0

    for X, Y in dataloader_test:
        X = X.reshape(-1, INPUT_SIZE).to(device) # l'envoie à device en corrigeant
                                            # si nécessaire le format
        Y = Y.float().to(device)
        #print(Y)
        Y_pred = model(X.float())  #calcul de la prédiction sur le modèle entrainé
        #print(Y_pred)


        # dans mon cas je veux les outputs pour calculer brier dessus

        # value, index


        N_SAMPLES += X.shape[0]  # nb of examples in the current batch

        # combien de bonnes prédictions? sachant que dans la database on ne donne que des cas
        # ou aa_origine == aa_destination

        # value, index
        _, predictions = torch.max(Y_pred, 1)
        # print(predictions)

        _, real_values = torch.max(Y, 1)
        # print(real_values)

        N_CORRECT += (predictions == real_values).sum().item() # .item(): access 1 value tensor

    print("Example of a test:")
    print(f"--Prediction:\n {Y_pred}")
    print(f"--Expected values:\n {Y}")
    print("")
    print("Index of the Softmax MAX")
    print(f"--Prediction: {predictions}")
    print(f"--Expected values: {real_values}")

    print("")
    acc = 100.0 * N_CORRECT / N_SAMPLES
    print(f"ACCURACY: {acc} %")

"""
FFNN:
python3 ffnn_multi_layer.py max_neighbour.csv > max_neighbour_multi_layers_$$.txt 2>&1


OVERFITTING TEST + REMEDIATION: 
https://datahacker.rs/018-pytorch-popular-techniques-to-prevent-the-overfitting-in-a-neural-networks/

https://rubikscode.net/2021/08/02/pytorch-for-beginners-building-neural-networks/
"""
import argparse
import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import random_split
import torch.utils.data as data_utils
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib as plt

parser = argparse.ArgumentParser()
parser.add_argument("path_data", help="dataset global path", type=str)
args = parser.parse_args()

N_NEIGHBOUR = 3
N_CHARACTER = 9   # pour le test


#######################################

print("________________________________________________")
print("             LEARNING IDENTITY TEST             ")
print("      aa_destination == max(aa_contextual)      ")
print("________________________________________________")

# device config
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')




# hyper parameters
N_CLASSES = N_CHARACTER
INPUT_SIZE = N_CHARACTER*(N_NEIGHBOUR +1) # not sure
N_HIDDEN_LAYERS = 2
HIDDEN_SIZE = 5  # comment le choisir ?
N_EPOCHS = 10
BATCH_SIZE = 5
LEARNING_RATE = 0.001
TEST_SIZE = 0.2
VALIDATION_FRACTION = 0.2
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
# print(f"dataset.n_samples, dataset.n_features: {dataset.n_samples}, {dataset.n_features}")

X_train, X_test, Y_train, Y_test = train_test_split(dataset.X,
                                                    dataset.Y,
                                                    test_size=TEST_SIZE,
                                                    random_state=RANDOM_STATE)

# first_data = dataset[0] # premier ex récupérer grace à la méthode __getitem__
# X_first_data, Y_first_data = first_data
# print(f"nbre de features dans dataset: {dataset.n_features}")



################################################################
# DATALOAD

dataset_train = data_utils.TensorDataset(X_train, Y_train)
validation_size = int(VALIDATION_FRACTION*len(dataset_train))
print(f"validation_size: {validation_size}")
train_size = len(dataset_train) - validation_size
print(f"train_size: {train_size}")

train_data, val_data = random_split(dataset_train, [train_size, validation_size])

train_loader = DataLoader(train_data, BATCH_SIZE, shuffle=True, num_workers=4, pin_memory=True)
val_loader = DataLoader(val_data, BATCH_SIZE*2, num_workers=4, pin_memory=True) # POURQUOI *2?
                                                                                # pin_memory?
print(iter(val_loader).next())






# dataset_test = data_utils.TensorDataset(X_test, Y_test)
# dataloader_test = DataLoader(dataset=dataset_test, batch_size=BATCH_SIZE, shuffle=True, num_workers=2)

print("")
print("____________________")
print("  DATASET/DATALOAD  ")
print("____________________")
print("")
print("DATASET with 9 amino-acids encoded from 0 to 8")
print(f"N_NEIGHBOUR : {N_NEIGHBOUR}")
print(f"N_EXAMPLES (TOTAL): {len(dataset)}")
# print(f"Example of one X_train, Y_train: {first_data}")
print(f"N_EXAMPLES_TRAIN: {len(dataset_train)}")
# print(f"N_EXAMPLES_TEST: {len(dataset_test)}")
print("(aa_destination == max(aa_contextual)) == True  # for all the examples")

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
class FFNN(nn.Module):
    """Simple Feed Forward Neural Network with n hidden layers"""
    def __init__(self, input_size, num_hidden_layers, hidden_size, out_size, accuracy_function):
        super().__init__()
        self.accuracy_function = accuracy_function
        
        # Create first hidden layer
        self.input_layer = nn.Linear(input_size, hidden_size)
        
        # Create remaining hidden layers
        self.hidden_layers = nn.ModuleList()
        for i in range(0, num_hidden_layers):
            self.hidden_layers.append(nn.Linear(hidden_size, hidden_size))
        
        # Create output layer
        self.output_layer = nn.Linear(hidden_size, out_size)

        # # Softmax
        self.m = nn.Softmax(dim=1)
        
    def forward(self, x):
        # Utilize hidden layers and apply activation function
        out = self.input_layer(x)
        out = F.relu(out)
        
        for layer in self.hidden_layers:
            out = layer(out)
            out = F.relu(out)
        
        # Get predictions
        out = self.output_layer(out) # pas de relu en sortie?
        out = self.m(out)
        return out
    
    def training_step(self, batch):
        # Load batch
        X, Y = batch
        X = X.reshape(-1, INPUT_SIZE).to(device)
        Y = Y.to(device)
        # Generate predictions
        Y_pred = self(X.float())
        # print(f"Y_pred train: {Y_pred}")
        # print(f"Y val size train: {Y.size()}")
        # print(f"Y_pred val size train: {Y_pred.size()}")
        
        # Calculate loss
        loss = F.cross_entropy(Y_pred, Y.float())
        return loss
    
    def validation_step(self, batch):
        # Load batch
        X, Y = batch
        X = X.reshape(-1, INPUT_SIZE).to(device)
        Y = Y.to(device)
        # Generate predictions
        Y_pred = self(X.float())
        # print(f"Y val size val: {Y.size()}")
        # print(f"Y_pred val size val: {Y_pred.size()}")
        
        # Calculate loss
        loss = F.cross_entropy(Y_pred, Y.float())

        # Calculate accuracy
        acc = self.accuracy_function(Y_pred, Y.float())
        
        return {'val_loss': loss, 'val_acc': acc}
        
    def validation_epoch_end(self, outputs):
        batch_losses = [x['val_loss'] for x in outputs]
        
        # Combine losses and return mean value
        epoch_loss = torch.stack(batch_losses).mean()
        
        # Combine accuracies and return mean value
        batch_accs = [x['val_acc'] for x in outputs]
        epoch_acc = torch.stack(batch_accs).mean()
        return {'val_loss': epoch_loss.item(), 'val_acc': epoch_acc.item()}
    
    def epoch_end(self, epoch, result):
        print("Epoch: {} - Validation Loss: {:.4f}, Validation Accuracy: {:.4f}".format( \
						epoch, result['val_loss'], result['val_acc']))







class ModelTrainer():
    def fit(self, epochs, learning_rate, model, train_loader, val_loader, opt_func=torch.optim.SGD):
        history = []
        optimizer = opt_func(model.parameters(), learning_rate)

        for epoch in range(epochs):
            # Training 
            for batch in train_loader:
                loss = model.training_step(batch)
                loss.backward()
                optimizer.step()
                optimizer.zero_grad()

            # Validation
            result = self._evaluate(model, val_loader)
            model.epoch_end(epoch, result)
            history.append(result)
            
        return history

    def _evaluate(self, model, val_loader):
        outputs = [model.validation_step(batch) for batch in val_loader]
        return model.validation_epoch_end(outputs)




def accuracy(outputs, labels):
    _, preds = torch.max(outputs, dim=1)
    _, real_values = torch.max(labels, 1)

    return torch.tensor(torch.sum(preds == real_values).item() / len(preds))
    # return torch.tensor(torch.sum(preds == labels).item() / len(preds))


def plot_history(history):
    losses = [x['val_loss'] for x in history]
    plt.plot(losses, '-x')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')

    accuracies = [x['val_acc'] for x in history]
    plt.plot(accuracies, '-x')
    plt.xlabel('epoch')
    plt.ylabel('accuracy')
    plt.title('Loss and Accuracy')












## modèle
model = FFNN(INPUT_SIZE, N_HIDDEN_LAYERS, HIDDEN_SIZE, out_size=N_CLASSES, accuracy_function=accuracy).to(device)
print(model)

model_trainer = ModelTrainer()

training_history = []

training_history += model_trainer.fit(N_EPOCHS, LEARNING_RATE, model, train_loader, val_loader)

plot_history(training_history)

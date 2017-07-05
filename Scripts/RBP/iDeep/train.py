"""Train the 3 models for iDeep experiment

Author: Mohammadamin Barekatain
Affiliation: TUM
"""

from data import data
from model import model_ideep
from keras.callbacks import EarlyStopping
from sklearn.metrics import roc_auc_score
from data import dtdir
import os
import pandas as pd
import numpy as np


outdir = "../../../Output/RBP/iDeep/"

def train_ideep(protein, features, save_path, seed=None):

    train, valid, test = data(protein, features, seed=seed)
    model = model_ideep(train[0], features, seed=seed)

    list_train = [train[0][feature] for feature in features]
    list_valid = [valid[0][feature] for feature in features]
    list_test = [test[0][feature] for feature in features]

    if seed is not None:
        np.random.seed(seed)

    earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=0)
    print('model training')
    model.fit(list_train, train[1], batch_size=100, epochs=20, verbose=0, validation_data=(list_valid, valid[1]), callbacks=[earlystopper])


    predictions = model.predict_proba(list_test)
    auc = roc_auc_score(test[1], predictions)
    print("Test AUC: ", auc)
    dt = pd.DataFrame( {"y_true": test[1], "y_pred_proba": predictions[:,0]} )
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    dt.to_csv(os.path.join(save_path, protein + ".csv"))

    return auc


if __name__ == "__main__":

    seed = 1337


    print ("*** Original iDeep ***")
    sum = []
    features = ["X_RG", "X_CLIP", "X_RNA", "motif", "seq"]
    for protein in os.listdir(dtdir):
        if protein != '.DS_Store':
            print("========", protein, "=========")
            sum.append(train_ideep(protein, features, save_path=os.path.join(outdir, 'iDeep_original'), seed = seed))
    sum = np.array(sum)
    print ('Average auc:', sum.mean())
    print ('Std auc:', sum.std())


    print("*** iDeep + spline transformation of positional feature ***")
    sum = []
    features = ["X_RG", "X_CLIP", "X_RNA", "motif", "seq", "positions_gam"]
    for protein in os.listdir(dtdir):
        if protein != '.DS_Store':
            print("========", protein, "=========")
            sum.append(train_ideep(protein, features, save_path=os.path.join(outdir, 'iDeep_scaler_position_gam'), seed=seed))
    sum = np.array(sum)
    print ('Average auc:', sum.mean())
    print ('Std auc:', sum.std())


    print ("*** iDeep + piecewise linear transformation of positional feature ***")
    sum = []
    features = ["X_RG", "X_CLIP", "X_RNA", "motif", "seq", "positions_nat"]
    for protein in os.listdir(dtdir):
        if protein != '.DS_Store':
            print("========", protein, "=========")
            sum.append(train_ideep(protein, features, save_path=os.path.join(outdir, 'iDeep_scaler_position_relu'), seed=seed))
    sum = np.array(sum)
    print ('Average auc:', sum.mean())
    print ('Std auc:', sum.std())


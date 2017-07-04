from data import *
from model import *
from keras.callbacks import ModelCheckpoint, EarlyStopping
from sklearn.metrics import roc_auc_score
from data import dtdir
import os
import pandas as pd
import datetime
import numpy as np
import json


def train_model(protein, features, save_path, save_model=False):
    hyper_params = {
        'X_RG': {'num_hidden_1': 256, 'num_hidden_2': 256, 'rg': 0}, #1e-2
        'X_CLIP': {'num_hidden_1': (256*3), 'num_hidden_2': (256*3), 'rg': 0},
        'X_RNA': {'num_hidden_1': 128, 'num_hidden_2': 128, 'rg': 0},
        'motif': {'num_hidden_1': 128, 'num_hidden_2': 128, 'rg': 0},
        'seq': {'num_filters': 128, 'dropout': 0.25, 'kernel_size': 7, 'pool_size': 3},
        'positions_gam': {'num_filters': 6}
    }

    train, valid, test = data_2(protein, features)
    model = model_3(train[0], features, hyper_params)

    list_train = [train[0][feature] for feature in features]
    list_valid = [valid[0][feature] for feature in features]
    list_test = [test[0][feature] for feature in features]

    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    checkpointer = ModelCheckpoint(filepath=os.path.join(save_path, "weights.hdf5"), verbose=1, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=7, verbose=0)
    print('model training')
    model.fit(list_train, train[1], batch_size = 128, epochs=30, verbose=2, validation_data=(list_valid, valid[1]),
              callbacks=[checkpointer, earlystopper])
    model.load_weights(os.path.join(save_path, "weights.hdf5"))

    predictions = model.predict_proba(list_test)
    auc = roc_auc_score(test[1], predictions)
    print("Test AUC: ", auc)
    dt = pd.DataFrame({"y_true": test[1], "y_pred_proba": predictions[:, 0]})
    dt.to_csv(os.path.join(save_path, protein + ".csv"))
    if save_model:
        model.save(os.path.join(save_path, protein + "_model.h5"))

    with open(os.path.join(save_path, 'hyper_params.json'), 'w') as fp:
        json.dump(hyper_params, fp)

    return auc


def train_multi_class_weight(protein, features, weight_range, save_path, seed, loss='binary_crossentropy', save_model=False):

    train, valid, test = data(protein, features, seed=seed)

    list_train = [train[0][feature] for feature in features]
    list_valid = [valid[0][feature] for feature in features]
    list_test = [test[0][feature] for feature in features]

    y_train = train[1]
    y_valid = valid[1]
    # y_test = test[1]

    if loss in ['hinge', 'squared_hinge', 'categorical_hinge']:
        y_train[y_train==0] = -1
        y_valid[y_valid==0] = -1
    #     from keras.utils.np_utils import to_categorical
    #     y_train = to_categorical(y_train, num_classes=None)
    #     y_valid = to_categorical(y_valid, num_classes=None)
        # y_test = to_categorical(y_test, num_classes=None)

    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    auc = []
    for w in weight_range:

        if seed is not None:
            np.random.seed(seed)

        if loss == 'binary_crossentropy':
            class_weight = {0: 1, 1: w}
        elif loss in ['hinge', 'squared_hinge', 'categorical_hinge']:
            class_weight = {-1: 1, 1: w}

        checkpointer = ModelCheckpoint(filepath=os.path.join(save_path, "weights.hdf5"), verbose=1, save_best_only=True)
        earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=0)


        print('model training ', class_weight)
        model = model_ideep(train[0], features, loss=loss, seed=seed)
        model.fit(list_train, y_train, batch_size=100, epochs=50, verbose=2, validation_data=(list_valid, y_valid),
              callbacks=[earlystopper, checkpointer], class_weight=class_weight)
        model.load_weights(os.path.join(save_path, "weights.hdf5"))


        predictions = model.predict_proba(list_test)
        # if loss=='binary_crossentropy':
        auc.append(roc_auc_score(test[1], predictions))
        dt = pd.DataFrame({"y_true": test[1], "y_pred_proba": predictions[:, 0]})
        # elif loss in ['hinge', 'squared_hinge', 'categorical_hinge']:
        #     # auc.append(roc_auc_score(test[1], predictions[:, 1]))
        #     # dt = pd.DataFrame({"y_true": test[1], "y_pred_proba": predictions[:, 1]})

        print("Test AUC with w=%s: %f" % (str(w), auc[-1]))
        dt.to_csv(os.path.join(save_path, protein + ("_%s" % (str(w))) + ".csv"))

        if save_model:
            model.save(os.path.join(save_path, protein + ("_%s" % (str(w))) + "_model.h5"))

    return np.array(auc)

def train_ideep(protein, features, save_path = "../results/ideep", save_model=False):

    train, valid, test = data(protein, features)
    model = model_ideep(train[0], features)

    list_train = [train[0][feature] for feature in features]
    list_valid = [valid[0][feature] for feature in features]
    list_test = [test[0][feature] for feature in features]

    earlystopper = EarlyStopping(monitor='val_loss', patience=7, verbose=0)
    print('model training')
    
    model.fit(list_train, train[1], batch_size=100, epochs=20, verbose=2, validation_data=(list_valid, valid[1]), callbacks=[earlystopper])
    predictions = model.predict_proba(list_test)
    auc = roc_auc_score(test[1], predictions)
    print("Test AUC: ", auc)
    dt = pd.DataFrame( {"y_true": test[1], "y_pred_proba": predictions[:,0]} )
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    dt.to_csv(os.path.join(save_path, protein + ".csv"))
    if save_model:
        model.save(os.path.join(save_path, protein + "_model.h5"))
    return auc

# def train_ideep_merged(proteins, features, save_path="../results/ideep_merged"):
#
#     train, valid, test = data_merged(proteins, features)
#     model = model_ideep(train[0], features)
#
#     list_train = [train[0][feature] for feature in features]
#     list_valid = [valid[0][feature] for feature in features]
#     list_test = [test[0][feature] for feature in features]
#
#     earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=0)
#     print 'model training'
#
#     model.fit(list_train, train[1], batch_size=100, epochs=20, verbose=1, validation_data=(list_valid, valid[1]),
#               callbacks=[earlystopper])
#     predictions = model.predict_proba(list_test)
#     auc = roc_auc_score(test[1], predictions)
#     print "Test AUC: ", auc
#     dt = pd.DataFrame({"y_true": test[1], "y_pred_proba": predictions[:, 0]})
#     if not os.path.isdir(save_path):
#         os.makedirs(save_path)
#     dt.to_csv(os.path.join(save_path, "merged.csv"))
#     model.save(os.path.join(save_path, "model.h5"))
 # train_ideep_merged(os.listdir(dtdir)[1:], features)

if __name__ == "__main__":
    sum = []
    exp_name = 'Model_3'
    features = ["X_RG", "X_CLIP", "X_RNA", "motif", "seq", "positions_gam"]
    i = datetime.datetime.now()
    time = "%s.%s_%s.%s" % (i.month, i.day, i.hour, i.minute)
    for protein in reversed(os.listdir(dtdir)):
        if protein != '.DS_Store':
            print("========", protein, "=========")
            sum.append(train_model(protein, features, save_path = os.path.join("../results/", exp_name, time)))
            # sum.append(train_ideep(protein, features, save_path=os.path.join("../results/", exp_name, time)))
    sum = np.array(sum)
    print ('Average auc:', sum.mean())
    print ('Std auc:', sum.std())
    # seed = 1337
    # sum = []
    # exp_name = 'multi_class_weight'
    # features = ["X_RG", "X_CLIP", "X_RNA", "motif", "seq", "positions_gam"]
    # i = datetime.datetime.now()
    # time = "%s.%s_%s.%s" % (i.month, i.day, i.hour, i.minute)
    # weight_range =  np.logspace(-2, 2, num=5) # np.array([0.2, 0.5, 1, 1.5, 2])
    # for protein in reversed(os.listdir(dtdir)[:14]):
    #     if protein != '.DS_Store':
    #         print("========", protein, "=========")
    #         sum.append(train_multi_class_weight(protein, features, weight_range = weight_range,
    #                                 save_path=os.path.join("../results/", exp_name+"_"+str(seed)+"_softmax", time), loss='binary_crossentropy', seed=seed))
    # sum = np.array(sum)
    # print('Average auc:', sum.mean(axis=0))
    # print('Std auc:', sum.std(axis=0))
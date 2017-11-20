"""data() returning arrays necessary for training and testing

Author: Mohammadamin Barekatain
Affiliation: TUM
"""

# Parts of this script has been copied from https://github.com/xypan1232/iDeep

import os
import numpy as np
import gzip
import random
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Imputer
import pandas as pd
from concise.preprocessing import encodeSplines

dtdir = '../../../Data/RBP/iDeep/raw'

files = {"X_GO": "matrix_GeneOntology.tab.gz", "X_KMER": "matrix_RNAkmers.tab.gz", "X_RG": "matrix_RegionType.tab.gz",
         "X_CLIP": "matrix_Cobinding.tab.gz", "X_RNA": "matrix_RNAfold.tab.gz", "positions_nat": "positions.csv", "positions_gam": "positions.csv"}


def read_positions(positions_file):

    dt = pd.read_csv(positions_file)
    positions = dt.as_matrix(dt.columns[:8])
    return np.sign(positions) * np.log10(np.abs(positions) + 1)


def read_seq(seq_file):

    def get_RNA_seq_concolutional_array(seq, motif_len=4):

        seq = seq.replace('U', 'T')
        alpha = 'ACGT'
        row = (len(seq) + 2 * motif_len - 2)
        new_array = np.zeros((row, 4))
        for i in range(motif_len - 1):
            new_array[i] = np.array([0.25] * 4)
            new_array[row - 1 - i] = np.array([0.25] * 4)

        for i, val in enumerate(seq):
            i = i + motif_len - 1
            if val not in 'ACGT':
                new_array[i] = np.array([0.25] * 4)
                continue
            try:
                index = alpha.index(val)
                new_array[i][index] = 1
            except:
                pdb.set_trace()
        return new_array

    with gzip.open(seq_file, 'r') as fp:
        seq_list = []
        seq = ''
        for line in fp:
            if str(line)[2] == '>':
                if len(seq):
                    seq_array = get_RNA_seq_concolutional_array(seq)
                    seq_list.append(seq_array)

                seq = ''
            else:
                seq = seq + str(line)[2:-3]
        if len(seq):
            seq_array = get_RNA_seq_concolutional_array(seq)
            seq_list.append(seq_array)

    return np.array(seq_list)


def preprocess_data(X_train, X_valid, X_test, n_bases):

    for key in X_train:
        if key == "seq":
            continue
        scaler = MinMaxScaler()
        if key in ["X_RNA", "motif"]:
            scaler = StandardScaler()
        if key in ["positions_nat", "positions_gam"]:
            imp = Imputer(missing_values='NaN', strategy='median', axis=0)
            imp.fit(X_train[key])
            X_train[key] = imp.transform(X_train[key])
            X_valid[key] = imp.transform(X_valid[key])
            X_test[key] = imp.transform(X_test[key])

        scaler.fit(X_train[key])
        X_train[key] = scaler.transform(X_train[key])
        X_valid[key] = scaler.transform(X_valid[key])
        X_test[key] = scaler.transform(X_test[key])

        if key == "positions_nat":
            X_train[key] = np.expand_dims(X_train[key], axis=2)
            X_valid[key] = np.expand_dims(X_valid[key], axis=2)
            X_test[key] = np.expand_dims(X_test[key], axis=2)

        if key == "positions_gam":
            # NOTE: concise now implements a transoformer API - concise.preprocessing.EncodeSplines with methods:
            # .fit
            # .predict
            # and operates on multiple features in a single array. Please use that one instead.
            X_train[key] = encodeSplines(X_train[key], n_bases=n_bases, start=0, end=1)
            X_valid[key] = encodeSplines(X_valid[key], n_bases=n_bases, start=0, end=1)
            X_test[key] = encodeSplines(X_test[key], n_bases=n_bases, start=0, end=1)

def split_training_validation(X, y, validation_size=0.2, shuffle=False):
    """split sampels based on balnace classes"""

    num_samples = y.shape[0]
    classes = np.unique(y)
    training_mask = np.ones(num_samples, dtype='bool')

    for cl in classes:
        cl_indices = np.where(y == cl)[0]
        if shuffle:
            random.shuffle(cl_indices)
        valid_cl_length = int(len(cl_indices) * validation_size)
        valid_cl_indices = cl_indices[:valid_cl_length]
        training_mask[valid_cl_indices] = False

    training_indices = np.where(training_mask == True)[0]
    validation_indices = np.where(training_mask == False)[0]
    random.shuffle(training_indices)
    random.shuffle(validation_indices)
    y_train = y[training_indices]
    y_valid = y[validation_indices]

    X_train, X_valid = ({}, {})
    for key in X:
        X_train[key] = X[key][training_indices]
        X_valid[key] = X[key][validation_indices]

    return X_train, y_train, X_valid, y_valid


def data(protein, features, n_bases=32, seed=None):
    if seed is not None:
        np.random.seed(seed)

    path = os.path.join(dtdir, protein, '5000')
    X, X_test = ({}, {})
    y = np.loadtxt(gzip.open(os.path.join(path, 'training_sample_0', "matrix_Response.tab.gz")), skiprows=1)
    y_test = np.loadtxt(gzip.open(os.path.join(path, 'test_sample_0', "matrix_Response.tab.gz")), skiprows=1)

    for feature in features:
        print('loading ', feature)
        if feature == "motif":
            X["motif"] = np.loadtxt(gzip.open(os.path.join(path, 'training_sample_0', 'motif_fea.gz')), skiprows=1, usecols=list(range(1, 103)))
            X_test["motif"] = np.loadtxt(gzip.open(os.path.join(path, 'test_sample_0', 'motif_fea.gz')), skiprows=1, usecols=list(range(1, 103)))
        elif feature == "seq":
            X["seq"] = read_seq(os.path.join(path, 'training_sample_0', 'sequences.fa.gz'))
            X_test["seq"] = read_seq(os.path.join(path, 'test_sample_0', 'sequences.fa.gz'))
        elif feature in ["positions_nat", "positions_gam"]:
            X[feature] = read_positions(os.path.join(path, 'training_sample_0', files[feature]))
            X_test[feature] = read_positions(os.path.join(path, 'test_sample_0', files[feature]))
        else:
            X[feature] = np.loadtxt(gzip.open(os.path.join(path, 'training_sample_0', files[feature])), skiprows=1)
            X_test[feature] = np.loadtxt(gzip.open(os.path.join(path, 'test_sample_0', files[feature])), skiprows=1)

    X_train, y_train, X_valid, y_valid = split_training_validation(X, y)
    X.clear()
    preprocess_data(X_train, X_valid, X_test, n_bases)
    return (X_train, y_train), (X_valid, y_valid), (X_test, y_test)

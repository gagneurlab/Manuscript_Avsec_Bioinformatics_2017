"""model() returning the keras model

Author: Mohammadamin Barekatain
Affiliation: TUM
"""

# Parts of this script has been copied from https://github.com/xypan1232/iDeep

from keras.models import Sequential, model_from_config
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.normalization import BatchNormalization
from keras.layers.local import LocallyConnected1D
from keras.layers.convolutional import Convolution2D, MaxPooling2D, Conv1D
from keras.layers.pooling import MaxPooling1D
from keras.layers import Merge
import numpy as np

def model_ideep(train_data, features, seed=None):
    if seed is not None:
        np.random.seed(seed)

    def get_cnn_gam(train, num_filters=6):
        print('configure gam network for', train.shape)

        ext_n_bases = train.shape[2]
        model = Sequential()
        model.add(LocallyConnected1D(filters=num_filters, kernel_size=1, input_shape=(8, ext_n_bases), activation='relu'))
        model.add(Flatten())
        model.add(Dense(64, activation='relu'))

        return model

    def get_fnn_nat(train, ext_n_bases=32, num_filters=6):
        print('configure network for', train.shape)

        model = Sequential()
        model.add(LocallyConnected1D(filters=ext_n_bases, kernel_size=1, input_shape=(8, 1), activation='relu'))
        model.add(LocallyConnected1D(filters=num_filters, kernel_size=1, activation='relu'))
        model.add(Flatten())
        model.add(Dense(64, activation='relu'))
        return model

    def get_rnn_fea(train, num_hidden=128):
        print('configure network for', train.shape)

        model = Sequential()
        model.add(Dense(num_hidden, input_shape=(train.shape[1],), activation='relu'))
        model.add(BatchNormalization())
        model.add(Dropout(0.5))
        model.add(Dense(num_hidden, input_dim=num_hidden, activation='relu'))
        model.add(BatchNormalization())
        model.add(Dropout(0.5))
        return model


    def get_cnn_network():
        print('configure cnn network')

        nbfilter = 102
        model = Sequential()
        model.add(Conv1D(activation="relu", input_shape=(107, 4), padding="valid", strides=1, filters=nbfilter, kernel_size=7))
        model.add(MaxPooling1D(pool_size=3))
        model.add(Dropout(0.5))
        model.add(Flatten())
        model.add(Dense(nbfilter, activation='relu'))
        model.add(Dropout(0.25))

        return model

    fea_num_hidden = {"X_GO": 2048, "X_KMER": 2048, "X_RG": 256, "X_CLIP": (256*3), "X_RNA": 128, "motif": 128, "positions_nat": 64}
    net_list = []
    for feature in features:
        if feature == "seq":
            net_list.append(get_cnn_network())
        elif feature == "positions_nat":
            net_list.append(get_fnn_nat(train_data[feature]))
        elif feature == "positions_gam":
            net_list.append(get_cnn_gam(train_data[feature]))
        else:
            net_list.append(get_rnn_fea(train_data[feature], fea_num_hidden[feature]))

    model = Sequential()
    model.add(Merge(net_list, mode='concat'))
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='rmsprop')

    return model


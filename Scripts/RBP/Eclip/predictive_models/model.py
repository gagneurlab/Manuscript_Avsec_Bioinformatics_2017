"""RBP model
"""
import numpy as np
from keras.optimizers import Adam, SGD
import keras.initializers as ki
import concise.losses as closs
import keras.layers as kl
import keras.regularizers as kr
from keras.models import Sequential, Model
import concise.initializers as ci
import concise.regularizers as cr
import concise.metrics as cm
import concise.layers as cl
import concise.activations as ca
from concise.utils import model_data
import keras.backend as K
from concise.optimizers import AdamWithWeightnorm, data_based_init
from keras.layers.advanced_activations import PReLU
from copy import deepcopy

def get_pool(pooling):
    pooling = deepcopy(pooling)
    pool_type = pooling["pool_type"]
    del pooling["pool_type"]
    if pool_type == "max":
        return [kl.MaxPooling1D(pool_size=pooling["pool_size"]), kl.Flatten()]
    elif pool_type == "mean":
        return [kl.AveragePooling1D(pool_size=pooling["pool_size"]), kl.Flatten()]
    elif pool_type == "weight":
        return [cl.SmoothPositionWeight(**pooling), cl.GlobalSumPooling1D()]
    else:
        raise ValueError("")


def model_simple(train_data,
                 kernel_size=22,
                 filters=64,
                 pooling_final={"pool_type": "max",
                                "pool_size": 4},  # weight pooling also aplicable
                 activation="relu",
                 lr=0.001):
    seq_length = train_data[0]["seq"].shape[1]
    n_tasks = train_data[1].shape[1]

    # model
    m = Sequential([
        cl.ConvDNA(filters=filters,
                   kernel_size=kernel_size,
                   activation=activation,
                   seq_length=seq_length),
        kl.MaxPooling1D(4),
        kl.BatchNormalization(axis=1),
        kl.Conv1D(filters=filters * 2,  # 128
                  kernel_size=4,
                  activation=activation),
        *get_pool(pooling_final),
        kl.BatchNormalization(),
        kl.Dropout(0.5),
        kl.Dense(64, activation=activation),
        kl.Dropout(0.5),
        kl.Dense(n_tasks, activation="sigmoid"),
    ])
    m.compile(optimizer=Adam(lr=lr),
              loss="binary_crossentropy",
              metrics=["acc"],)
    return m


def pos_module(pos_length, ext_n_bases, ext_filters, feat_name, ext_pos_kwargs):
    """Used position module function, calls pos_spline_trans or pos_affine_relu
    """
    if ext_pos_kwargs["type"] == "gam":
        return pos_spline_trans(pos_length, ext_n_bases, ext_filters, feat_name, ext_pos_kwargs)
    elif ext_pos_kwargs["type"] == "relu":
        return pos_affine_relu(pos_length, ext_n_bases, ext_filters, feat_name, ext_pos_kwargs)
    else:
        raise ValueError("pos_type invalid")


def pos_spline_trans(pos_length, ext_n_bases, ext_filters, feat_name, ext_pos_kwargs):
    """Get the spline transformation module

    Returns the input and output node
    """
    reg = cr.GAMRegularizer(ext_n_bases,
                            l2_smooth=ext_pos_kwargs.get("l2_smooth", 0))
    inp = cl.InputSplines(pos_length, n_bases=ext_n_bases,
                          name="dist_" + feat_name)
    x = cl.ConvSplines(filters=ext_filters,
                       kernel_regularizer=reg,
                       name="conv_dist_" + feat_name)(inp)
    if not ext_pos_kwargs["as_track"]:
        x = kl.Flatten()(x)
    return inp, x


def pos_affine_relu(pos_length, ext_n_bases, ext_filters, feat_name, ext_pos_kwargs):
    """Get the affine+relu transformation module

    Returns the input and output node
    """
    inp = cl.Input((pos_length, 1),
                   name="raw_dist_" + feat_name)
    x = kl.Conv1D(filters=ext_n_bases,
                  kernel_size=1,
                  activation="relu")(inp)
    x = kl.Conv1D(filters=ext_filters, kernel_size=1,
                  name="conv_dist_" + feat_name)(x)
    if not ext_pos_kwargs["as_track"]:
        x = kl.Flatten()(x)
    return inp, x


# TODO - add secondary structure
# structure argument
def model(train_data,
          activation="relu",
          kernel_size=10,
          filters=16,
          conv2_use_skip=False,
          internal_pos={"name": "global_maxpool"},
          # {"name": "strided_maxpool", "pool_size": 3}
          # {"name": "maxpool+rnn_sequence", "pool_size": 3, "dropout": 0.1}
          # {"name": "rnn", "pool_size": 3, "dropout": 0.1}
          # {"name": "maxpool+weight_sum", "pool_size": 3, "n_bases": 10, "share_splines": False,
          # "l2_smooth": 1e-5, "l2": 0}
          external_pos=None,  # None, {"type": "gam", "as_track": True, "units": 1}
          dropout_rate=0.5,
          n_hidden=100,
          use_batchnorm=False,
          use_weightnorm=False,
          lr=0.001):
    """Returns keras model for modelling arrays returned by data.data()
    """
    # config
    seq_length = train_data[0]["seq"].shape[1]
    n_tasks = train_data[1].shape[1]
    ext_n_bases = train_data[0]["dist_tss"].shape[2]
    activation = PReLU() if activation == "PReLU" else activation

    inputs = []
    # position module
    # ---------------
    if external_pos is not None:
        # conf
        external_pos["as_track"] = train_data[0]["dist_tss"].shape[1] != 1
        if external_pos["as_track"]:
            pos_length = seq_length - kernel_size + 1
        else:
            pos_length = 1

        ext_filters = external_pos["units"]

        pos_inputs, pos_outputs = tuple(zip(*[pos_module(pos_length=pos_length,
                                                         ext_n_bases=ext_n_bases,
                                                         ext_filters=ext_filters,
                                                         feat_name=feat_name,
                                                         ext_pos_kwargs=external_pos)
                                              for feat_name in external_pos.get("feat_names",
                                                                                ["gene_start", "gene_end"])]))
        inputs += list(pos_inputs)
        pos_outputs = list(pos_outputs)

    # sequence module
    # ----------------
    # initialize conv kernels to known motif pwm's
    seq_input = cl.InputDNA(seq_length, name="seq")
    inputs += [seq_input]
    x_pwm = cl.ConvDNA(filters=filters,
                       kernel_size=kernel_size,
                       activation=activation
                       )(seq_input)
    if use_batchnorm:
        x_pwm = kl.BatchNormalization(axis=1)(x_pwm)

    # inject external_pos as a track
    if external_pos is not None and external_pos["as_track"]:
        x_pwm = kl.concatenate([x_pwm] + pos_outputs, axis=-1)

    x = kl.Conv1D(filters, kernel_size=1,
                  activation=activation
                  )(x_pwm)
    if conv2_use_skip:
        x = kl.concatenate([x_pwm, x])  # skip connection ?
    if use_batchnorm:
        x = kl.BatchNormalization(axis=1)(x)

    # summarize across sequence
    # -------------------------
    if internal_pos["name"] == "global_maxpool":
        x = kl.GlobalMaxPool1D()(x)
    elif internal_pos["name"] == "strided_maxpool":
        x = kl.MaxPool1D(pool_size=internal_pos["pool_size"])(x)
        x = kl.Flatten()(x)
    elif internal_pos["name"] == "maxpool+rnn_sequence":
        x = kl.MaxPool1D(pool_size=internal_pos["pool_size"])(x)
        x = kl.Bidirectional(kl.LSTM(filters,
                                     dropout_W=internal_pos["dropout"],
                                     dropout_U=internal_pos["dropout"],
                                     return_sequences=True)
                             )(x)
        x = kl.Flatten()(x)
    elif internal_pos["name"] == "rnn":
        x = kl.MaxPool1D(pool_size=internal_pos["pool_size"])(x)
        x = kl.Bidirectional(kl.LSTM(filters,
                                     dropout_W=internal_pos["dropout"],
                                     dropout_U=internal_pos["dropout"])
                             )(x)
    elif internal_pos["name"] == "maxpool+weight_sum":
        x = kl.MaxPool1D(pool_size=internal_pos["pool_size"])(x)
        x = cl.SmoothPositionWeight(n_bases=internal_pos.get("n_bases", 10),
                                    share_splines=internal_pos.get("share_splines", False),
                                    l2_smooth=internal_pos.get("l2_smooth", 0),
                                    l2=internal_pos.get("l2", 0)
                                    )(x)
        x = cl.GlobalSumPooling1D()(x)
    else:
        raise ValueError("invalid internal_pos")
    if use_batchnorm:
        x = kl.BatchNormalization()(x)
    x = kl.Dropout(dropout_rate)(x)

    # append external_pos as a scalar
    # -------------------------------
    if external_pos is not None and not external_pos["as_track"]:
        x = kl.concatenate([x] + pos_outputs, axis=-1)

    # FC layers
    # ---------
    x = kl.Dense(n_hidden, activation=activation)(x)
    if use_batchnorm:
        x = kl.BatchNormalization()(x)
    x = kl.Dropout(dropout_rate)(x)
    outputs = kl.Dense(n_tasks, activation="sigmoid")(x)

    # compile model
    # -------------
    m = Model(inputs, outputs)
    if use_weightnorm:
        optimizer = AdamWithWeightnorm(lr=lr)
    else:
        optimizer = Adam(lr=lr)

    m.compile(optimizer=optimizer, loss="binary_crossentropy",
              metrics=["acc"])

    if use_weightnorm:
        data_based_init(m, model_data.subset(train_data, np.arange(500))[0])
    return m

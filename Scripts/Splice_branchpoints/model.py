"""Branchpoint model
"""
import numpy as np
from keras.optimizers import Adam
import keras.initializers as ki
import concise.losses as closs
import keras.layers as kl
import keras.regularizers as kr
from keras.engine import Model
import concise.initializers as ci
import concise.regularizers as cr
import concise.metrics as cm
import concise.layers as cl
import concise.activations as ca
from concise.utils import model_data
import keras.backend as K
from concise.optimizers import AdamWithWeightnorm, data_based_init
from keras.layers.advanced_activations import PReLU


def model(train_data,
          nonlinearity="relu",  # 'relu' or 'exp' - TODO - PRelu?
          filters=1,
          init_motifs={"use_pssm": True,
                       "stddev": 0.1,
                       "mean_max_scale": 0.0,
                       },
          # init_sd_w=1e-3,
          pos_effect={"type": "gam",  # alternative relu - TODO - add activation for relu
                      "l2_smooth": 1e-5,
                      "l2": 1e-5,
                      "activation": None,
                      "use_bias": False,
                      "merge": {"type": "multiply"},
                      },
          use_weightnorm=False,
          l1_weights=0,
          l1_motif=0,
          hidden_fc=None,  # {"n_hidden": 4, activation="PReLU", "l1": 1e-4, "l2": 1e-5, "dropout_rate": 0.5}
          # learning rate
          lr=0.002
          ):

    seq_length = train_data[0]["seq"].shape[1]
    n_tasks = 1 if len(train_data[1].shape) == 1 else train_data[1].shape[1]
    assert seq_length - (11 - 1) == n_tasks
    pwm_list = train_data[4]
    pos_features = train_data[3]
    kernel_size = 11  # hard-coded for comparison with branchpointer

    # sequence module
    # ----------------
    # initialize conv kernels to known motif pwm's
    if init_motifs is not None:
        if init_motifs.get("n_pwm", None) is not None:
            pwm_list = pwm_list[:init_motifs["n_pwm"]]
        if init_motifs["use_pssm"]:
            kernel_initializer = ci.PSSMKernelInitializer(pwm_list, stddev=init_motifs["stddev"])
            bias_initializer = ci.PSSMBiasInitializer(pwm_list, kernel_size=kernel_size,
                                                      mean_max_scale=init_motifs.get("mean_max_scale", 0.0))
        else:
            kernel_initializer = ci.PWMKernelInitializer(pwm_list, stddev=init_motifs["stddev"])
            bias_initializer = ci.PWMBiasInitializer(pwm_list, kernel_size=kernel_size,
                                                     mean_max_scale=init_motifs.get("mean_max_scale", 0.0))
    else:
        kernel_initializer = "glorot_uniform"
        bias_initializer = "zeros"

    seq_input = cl.InputDNA(seq_length, name="seq")
    # convolution
    activation = PReLU() if nonlinearity == "PReLU" else nonlinearity

    x = cl.ConvDNA(filters=filters,
                   kernel_size=kernel_size,
                   kernel_regularizer=kr.l1(l1_motif),  # Regularization
                   activation=activation,
                   kernel_initializer=kernel_initializer,
                   bias_initializer=bias_initializer
                   )(seq_input)

    # positional module
    # -----------------

    # optional positional effect
    if pos_effect is not None:
        # config
        if pos_effect["merge"]["type"] == "multiply":
            pos_filters = filters
            merge_fun = kl.multiply
        elif pos_effect["merge"]["type"] == "concatenate":
            # if we concatenate, then the number of filters should be 1
            pos_filters = 1
            merge_fun = kl.concatenate
        elif pos_effect["merge"]["type"] == "add":
            pos_filters = filters
            merge_fun = kl.add
        else:
            raise ValueError("pos_effect[\"merge\"][\"type\"] needs to be from {multiply, concatenate, add}")

        pos_in_layers = []
        pos_out_layers = []
        for pos_feat in pos_features:
            pos_effect_type = pos_effect.get("type", "gam")
            if pos_effect_type == "gam":
                pos_effect["n_bases"] = train_data[0]["dist2"].shape[2]
                pos_in = cl.InputDNAQuantitySplines(n_tasks,  # valid padding from conv
                                                    n_bases=pos_effect["n_bases"],
                                                    name=pos_feat)
                pos_out = cl.ConvDNAQuantitySplines(filters=pos_filters,
                                                    kernel_regularizer=cr.GAMRegularizer(n_bases=pos_effect["n_bases"],
                                                                                         l2_smooth=pos_effect["l2_smooth"],
                                                                                         l2=pos_effect["l2"]),
                                                    use_bias=pos_effect["use_bias"],
                                                    activation=pos_effect.get("activation", None),
                                                    name=pos_feat + "_conv"
                                                    )(pos_in)
            elif pos_effect_type == "relu":
                assert train_data[0]["dist2"].shape[2] == 1
                pos_in = kl.Input((n_tasks, 1), name=pos_feat)
                pos_features = kl.Conv1D(filters=pos_effect["n_bases"],
                                         kernel_size=1,
                                         use_bias=True,
                                         activation="relu",
                                         name=pos_feat + "_conv_features"
                                         )(pos_in)
                pos_out = kl.Conv1D(filters=pos_filters,
                                    kernel_size=1,
                                    kernel_regularizer=kr.L1L2(l2=pos_effect["l2"]),
                                    use_bias=pos_effect["use_bias"],
                                    activation=pos_effect.get("activation", None),
                                    name=pos_feat + "_conv_combine"
                                    )(pos_features)
            else:
                raise ValueError("pos_effect[\"type\"] is invalid")
            if pos_effect["merge"]["type"] == "multiply":
                pos_out = kl.Lambda(lambda x: 1.0 + x)(pos_out)
            # TODO - implement splines sharing?
            pos_in_layers.append(pos_in)
            pos_out_layers.append(pos_out)

        # merge the layers
        # ----------------
        x = merge_fun([x] + pos_out_layers)
        input_list = [seq_input] + pos_in_layers
        if pos_effect["merge"]["type"] == "concatenate" and \
           pos_effect["merge"].get("hidden_fc", None) is not None:
            hidden_fc = pos_effect["merge"]["hidden_fc"]
            act_string = hidden_fc.get("activation", "relu")
            activation = PReLU() if act_string == "PReLU" else act_string
            for i in range(hidden_fc["n_layers"]):
                x = kl.Conv1D(hidden_fc["n_hidden"],
                              # TODO - maybe use more of a segmentatino taste with kernel_size > 1
                              kernel_size=1,
                              activation="relu",
                              kernel_regularizer=kr.L1L2(l1=hidden_fc.get("l1", 0), l2=hidden_fc.get("l2", 0)))(x)
                if hidden_fc["dropout_rate"]:  # non-zero
                    x = kl.Dropout(hidden_fc["dropout_rate"])(x)
    else:
        input_list = seq_input

    if pos_effect is None and filters == 1:
        predictions = kl.Activation("sigmoid")(x)
    else:
        predictions = kl.Conv1D(filters=1,
                                kernel_size=1,
                                kernel_regularizer=kr.l1(l1_weights),
                                activation="sigmoid",
                                name="Conv1D_final_layer"
                                # kernel_initializer=ki.RandomNormal(stddev=init_sd_w) # TODO - which initialization?
                                )(x)

    # predictions = kl.Flatten()(predictions)  # remove the last dimention
    model = Model(input_list, predictions)

    if use_weightnorm:
        optimizer = AdamWithWeightnorm(lr=lr)
    else:
        optimizer = Adam(lr=lr)

    model.compile(optimizer=optimizer, loss=closs.binary_crossentropy_masked,
                  metrics=[cm.accuracy, cm.f1,
                           cm.precision, cm.recall, cm.sensitivity,
                           cm.specificity, cm.fdr],
                  sample_weight_mode="temporal")

    if use_weightnorm:
        data_based_init(model, model_data.subset(train_data, np.arange(500))[0])
    return model

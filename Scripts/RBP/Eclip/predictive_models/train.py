#!/usr/bin/env python
"""Train the models
"""
from hyperopt import fmin, tpe, hp, pyll
from concise.hyopt import CompileFN, CMongoTrials, test_fn
from copy import deepcopy
import numpy as np
from glob import glob
import os
import data
import model

PROC_DIR = "/s/project/deepcis/encode/eclip/processed"
MAX_EVALS = 50

KILL_TIMEOUT = 60 * 20  # 20 minutes

RBP_LIST = ["UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2"]

# all rbp's
RBP_ALL = [os.path.basename(x).replace(".csv", "")
           for x in glob(PROC_DIR + "/design_matrix/train/*.csv")]

DB_NAME = "RBP__Eclip"

RUN_TEST = False

# functions
# ---------


def print_exp(exp_name):
    print("-" * 40 + "\nexp_name: " + exp_name)


def eval_expr(expr):
    a = np.array([pyll.stochastic.sample(expr) for x in range(100)])
    return np.unique(a)


def run_DeepNN_trials(exp_name, fn, hyper_params,
                      run_test=True,
                      max_evals=MAX_EVALS, rbp_list=RBP_LIST):
    """run DeepNN trials with custom pos_as_track and external_pos arguments
    """
    print_exp(exp_name)
    # -----
    if run_test:
        # Test all of them
        for rbp in rbp_list[:1]:
            # config
            c_hyper_params = deepcopy(hyper_params)
            c_hyper_params["data"]["rbp_name"] = rbp
            c_exp_name = exp_name + "_" + rbp
            fn.exp_name = c_exp_name
            # run
            test_fn(fn, c_hyper_params)

    for rbp in rbp_list:
        # config
        c_hyper_params = deepcopy(hyper_params)
        c_hyper_params["data"]["rbp_name"] = rbp
        c_exp_name = exp_name + "_" + rbp
        fn.exp_name = c_exp_name
        print_exp(c_exp_name)
        # run
        trials = CMongoTrials(DB_NAME, c_exp_name, kill_timeout=30 * 60)
        best = fmin(fn, c_hyper_params, trials=trials, algo=tpe.suggest, max_evals=max_evals)
        print("best_parameters: " + str(best))


if __name__ == "__main__":

    # --------------------------------------------
    exp_name = "DeepNN"
    print_exp(exp_name)
    # -----
    fn = CompileFN(DB_NAME, exp_name,
                   data_fn=data.data,
                   model_fn=model.model,
                   add_eval_metrics=["auprc", "auc"],
                   loss_metric="auprc",
                   loss_metric_mode="max",
                   valid_split=None)  # use it from the data.data

    hyper_params = {
        "data": {"n_bases": hp.choice("d_n_bases", (10, 20, 30)),
                 "pos_class_weight": 1.0,
                 "pos_as_track": False,
                 },
        "shared": {"kernel_size": 10},
        "model": {"activation": "relu",
                  "filters": 16,
                  "internal_pos": hp.choice("m_internal_pos", (
                      {"name": "global_maxpool"},
                      {"name": "strided_maxpool", "pool_size": 4},
                      {"name": "maxpool+weight_sum", "pool_size": 4,
                       "n_bases": 10, "share_splines": False,
                       "l2_smooth": hp.loguniform("m_l2_smooth",
                                                  np.log(1e-12),
                                                  np.log(1e-1)),
                       "l2": 0}
                  )),
                  "external_pos": None,
                  "dropout_rate": hp.uniform("m_dropout_rate", 0, 0.7),
                  "use_batchnorm": hp.choice("m_use_batchnorm", (False, True)),
                  "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(1e-2)),
                  },
        "fit": {"epochs": 150,
                "patience": 5,
                "batch_size": 128,
                "use_weight": False
                }
    }

    c_hyper_params = deepcopy(hyper_params)
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=RUN_TEST,
                      max_evals=MAX_EVALS,
                      rbp_list=RBP_LIST)

    # --------------------------------------------
    exp_name = "DeepNN_scalar_position_gam"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["data"]["pos_as_track"] = False
    c_hyper_params["model"]["external_pos"] = {"type": "gam",
                                               "scale": "log",
                                               "units": 1}

    # run for all the RBP's
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=RUN_TEST,
                      max_evals=MAX_EVALS,
                      rbp_list=RBP_LIST)
    # --------------------------------------------
    exp_name = "DeepNN_scalar_position_relu"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["data"]["pos_as_track"] = False
    c_hyper_params["model"]["external_pos"] = {"type": "relu",
                                               "scale": "log",
                                               "units": 1}
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=RUN_TEST,
                      max_evals=MAX_EVALS,
                      rbp_list=RBP_LIST)

    # --------------------------------------------
    exp_name = "DeepNN_track_position_gam"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["data"]["pos_as_track"] = True
    c_hyper_params["model"]["external_pos"] = {"type": "gam",
                                               "scale": "log",
                                               "units": 1}
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=RUN_TEST,
                      max_evals=MAX_EVALS,
                      rbp_list=RBP_LIST)

    # --------------------------------------------
    exp_name = "DeepNN_track_position_relu"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["data"]["pos_as_track"] = True
    c_hyper_params["model"]["external_pos"] = {"type": "relu",
                                               "scale": "log",
                                               "units": 1}

    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=RUN_TEST,
                      max_evals=MAX_EVALS,
                      rbp_list=RBP_LIST)

    # --------------------------------------------
    exp_name = "DeepNN_scalar_position_gam_all"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["data"]["pos_as_track"] = False
    c_hyper_params["model"]["external_pos"] = {"type": "gam",
                                               "scale": "log",
                                               "units": 1}
    c_hyper_params["model"]["kernel_size"] = hp.choice("m_kernel_size",
                                                       (10, 15))
    c_hyper_params["model"]["filters"] = hp.choice("m_filter",
                                                   (16, 32, 64, 128))

    # run for all the RBP's
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=RUN_TEST,
                      max_evals=20,
                      rbp_list=RBP_ALL)

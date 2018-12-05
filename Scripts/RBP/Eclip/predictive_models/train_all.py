#!/usr/bin/env python
"""Train and tune hyper-parameters for all 112 RBPs and 8 positional features

Author: Ziga Avsec
Affiliation: TUM
"""
from hyperopt import fmin, tpe, hp, pyll
from copy import deepcopy
import numpy as np
from glob import glob
import os
import data
import model
from joblib import Parallel, delayed
from mongodb_setup import host, port
import argparse
from concise.hyopt import CompileFN, CMongoTrials, test_fn


def print_exp(exp_name):
    print("-" * 40 + "\nexp_name: " + exp_name)


DIR_ROOT = data.DIR_ROOT
PROC_DIR = DIR_ROOT + "/processed"

PROC_DIR = DIR_ROOT + "/processed"
MAX_EVALS = 20
MAX_EVALS_subset = 50

KILL_TIMEOUT = 60 * 80  # 30 minutes

RBP_LIST = ["UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2"]

# all rbp's
RBP_ALL = [os.path.basename(x).replace(".csv", "")
           for x in glob(PROC_DIR + "/design_matrix/train/*.csv") if "extended" not in x]

RUN_TEST = False

DB_NAME = "RBP__Eclip"
POS_FEATURES = ['tss', 'polya', 'exon_intron', 'intron_exon', 'start_codon',
                'stop_codon', 'gene_start', 'gene_end']

POS_FEATURES_start_end = ["gene_start", "gene_end"]


# functions
# ---------
class RunFN():

    def __init__(self, exp_name, fn, hyper_params, max_evals, db_name=DB_NAME):
        self.exp_name = exp_name
        self.fn = fn
        self.hyper_params = hyper_params
        self.max_evals = max_evals
        self.db_name = db_name

    def __call__(self, rbp):
        # config
        c_hyper_params = deepcopy(self.hyper_params)
        c_hyper_params["data"]["rbp_name"] = rbp
        c_exp_name = exp_name + "_" + rbp
        fn.exp_name = c_exp_name
        print_exp(c_exp_name)
        # run
        trials = CMongoTrials(self.db_name, c_exp_name,  ip=host,
                              kill_timeout=KILL_TIMEOUT, port=port)
        best = fmin(fn, c_hyper_params, trials=trials,
                    algo=tpe.suggest, max_evals=self.max_evals)
        print("best_parameters: " + str(best))


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

    # Run in parallel
    results = Parallel(n_jobs=-1, backend="threading")(
        map(delayed(RunFN(exp_name, fn, hyper_params, max_evals)), rbp_list))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="run hyper-parameter optimization")
    parser.add_argument('--notest', action='store_true')
    args = parser.parse_args()

    run_test = not args.notest
    rbp_set = RBP_ALL
    rbp_subset = ["UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2"]

    print(rbp_set)
    # --------------------------------------------
    exp_name = "DeepNN_ext"
    print_exp(exp_name)
    # -----
    fn = CompileFN(DB_NAME, exp_name,
                   data_fn=data.data_extended_cached,
                   model_fn=model.model,
                   add_eval_metrics=["auprc", "auc"],
                   loss_metric="auprc",
                   loss_metric_mode="max",
                   save_dir="hyperopt_train_all",
                   valid_split=None)  # use it from the data.data

    hyper_params = {
        "data": {"n_bases": hp.choice("d_n_bases", (10, 20, 30)),
                 "pos_class_weight": 1.0,
                 "valid_chr": [1, 3],
                 "test_chr": [2, 4, 6, 8, 10]
                 },
        "model": {"activation": "relu",
                  "filters": 16,
                  "kernel_size": 11,
                  "internal_pos": {"name": "strided_maxpool", "pool_size": 4},
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
    chp_test = deepcopy(hyper_params)
    chp_test["data"]["rbp_name"] = "UPF1"
    # test_fn(fn, chp_test)
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=run_test,
                      max_evals=MAX_EVALS,
                      rbp_list=rbp_set)

    # --------------------------------------------
    exp_name = "DeepNN"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["model"]["external_pos"] = None

    # run for all the RBP's
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=run_test,
                      max_evals=MAX_EVALS,
                      rbp_list=rbp_subset)
    # --------------------------------------------
    exp_name = "DeepNN_scalar_position_gam"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["model"]["external_pos"] = {"type": "gam",
                                               "feat_names": POS_FEATURES_start_end,
                                               "units": 1}

    # run for all the RBP's
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=run_test,
                      max_evals=MAX_EVALS,
                      rbp_list=rbp_subset)
    # --------------------------------------------
    exp_name = "DeepNN_scalar_position_ext_gam"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["model"]["external_pos"] = {"type": "gam",
                                               "feat_names": POS_FEATURES,
                                               "units": 1}

    # run for all the RBP's
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=run_test,
                      max_evals=MAX_EVALS,
                      rbp_list=rbp_set)
    # --------------------------------------------
    exp_name = "DeepNN_scalar_position_ext_relu"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["model"]["external_pos"] = {"type": "relu",
                                               "feat_names": POS_FEATURES,
                                               "units": 1}
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=run_test,
                      max_evals=MAX_EVALS,
                      rbp_list=rbp_set)
    # --------------------------------------------
    exp_name = "DeepNN_track_position_ext_gam"
    # -----

    c_hyper_params = deepcopy(hyper_params)
    c_hyper_params["data"]["pos_as_track"] = True
    c_hyper_params["model"]["external_pos"] = {"type": "gam",
                                               "feat_names": POS_FEATURES,
                                               "units": 1}
    run_DeepNN_trials(exp_name, fn, c_hyper_params,
                      run_test=run_test,
                      max_evals=MAX_EVALS,
                      rbp_list=rbp_set)

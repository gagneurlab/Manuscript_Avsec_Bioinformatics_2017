#!/usr/bin/env python
from concise.hyopt import CMongoTrials, get_data
import numpy as np
import pandas as pd
import data
from glob import glob
import os

from joblib import Parallel, delayed

def print_exp(exp_name):
    print("-" * 40 + "\nexp_name: " + exp_name)


def get_predictions_method(exp_name, rbp):
    cache = False
    data_fn = data.data_extended
    # -----
    c_exp_name = exp_name + "_" + rbp
    print_exp(c_exp_name)

    basedir = "{root}/processed/predictions/{rbp}/".format(root=DIR_ROOT, rbp=rbp)
    out_csv = basedir + "{method}.csv".format(method=exp_name)
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    if cache and os.path.exists(out_csv):
        return
    # get trials
    trials = CMongoTrials(DB_NAME, c_exp_name, ip=HOST, kill_timeout=30 * 60)

    # get best trial parameters
    tid = trials.best_trial_tid()
    param = trials.get_param(tid)

    # load_model
    m = trials.load_model(tid)

    # load_data
    train, valid, test = get_data(data_fn, param)

    # predict
    y_pred = m.predict(test[0], verbose=0)
    y_true = test[1]

    # save
    dt_pred = pd.DataFrame({"y_true": y_true.reshape((-1,)),
                            "y_pred": y_pred.reshape((-1,))})

    dt_pred.to_csv(out_csv)


if __name__ == "__main__":

    PROC_DIR = "/s/project/deepcis/encode/eclip/processed"

    DIR_ROOT = "/s/project/deepcis/encode/eclip/"
    RBP_LIST = ["UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2"]

    # all rbp's
    RBP_ALL = [os.path.basename(x).replace(".csv", "")
               for x in glob(PROC_DIR + "/design_matrix/train/*.csv") if "_extended" not in x]

    DB_NAME = "RBP__Eclip"
    HOST = "ouga03"

    # EXPERIMENTS = ["DeepNN",
    #                "DeepNN_scalar_position_gam",
    #                "DeepNN_scalar_position_relu",
    #                "DeepNN_track_position_gam",
    #                "DeepNN_track_position_relu"]
    EXPERIMENTS = ["DeepNN_ext",
                   "DeepNN_scalar_position_ext_gam",
                   "DeepNN_scalar_position_ext_relu"]

    with Parallel(n_jobs=30) as parallel:
        parallel(delayed(get_predictions_method)(exp_name, rbp)
                 for exp_name in EXPERIMENTS
                 for rbp in RBP_ALL)

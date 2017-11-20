#!/usr/bin/env python
"""Retrieve the test-set predictions for the best model from each experiment

Author: Ziga Avsec
Affiliation: TUM
"""
import numpy as np
import pandas as pd
import data
from glob import glob
from concise.hyopt import CMongoTrials, get_data
import os
from train_all import DIR_ROOT, PROC_DIR, RBP_LIST, RBP_ALL, RBP_rerun, DB_NAME, HOST
import logging
import argparse
from pprint import pprint
logging.getLogger().setLevel(logging.WARN)

from joblib import Memory, Parallel, delayed

memory = Memory(cachedir=data.CACHE_DIR, verbose=0)


def print_exp(exp_name, rbp):
    print("-" * 40 + "\nexp: {0}; rbp: {1}".format(exp_name, rbp))


def get_trials(exp_name, rbp):
    return CMongoTrials(DB_NAME, exp_name + "_" + rbp, ip=HOST)


def get_predictions_overall(exp_name, rbp, ignore_cache=False):
    # the following function is cached if for the same experiment,
    # the same tid is the best one
    save_predictions_cached(exp_name, rbp,
                            best_tid=get_trials(exp_name, rbp).best_trial_tid(),
                            ignore_cache=ignore_cache)


@memory.cache(ignore=['ignore_cache'])
def save_predictions_cached(exp_name, rbp, best_tid, ignore_cache=False):
    data_fn = data.data_extended_cached
    print_exp(exp_name, rbp)

    basedir = "{root}/processed/predictions/{rbp}/".format(root=DIR_ROOT, rbp=rbp)
    out_csv = basedir + "{method}.csv".format(method=exp_name)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    trials = get_trials(exp_name, rbp)

    # no trials yet - return None
    if trials.n_ok() == 0:
        return

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
    return out_csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run predictions")
    parser.add_argument("--ignore_cache", action="store_true",
                        help="If enabled, ignore the cached computation")
    parser.add_argument("--n_jobs", type=int, default=10,
                        help="number of jobs to launch in parallel")
    args = parser.parse_args()

    EXPERIMENTS = {"DeepNN_ext": RBP_ALL,
                   "DeepNN": RBP_ALL,
                   "DeepNN_scalar_position_gam": RBP_ALL,
                   "DeepNN_scalar_position_ext_gam": RBP_ALL,
                   "DeepNN_scalar_position_ext_relu": RBP_ALL,
                   "DeepNN_track_position_ext_gam": RBP_ALL}

    pprint("Running for experiments:")
    pprint(list(EXPERIMENTS))
    with Parallel(n_jobs=args.n_jobs) as parallel:
        parallel(delayed(get_predictions_overall)(exp_name, rbp,
                                                  ignore_cache=args.ignore_cache)
                 for exp_name, rbps in EXPERIMENTS.items()
                 for rbp in rbps)

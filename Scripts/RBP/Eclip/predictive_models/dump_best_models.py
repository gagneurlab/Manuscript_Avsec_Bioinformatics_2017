#!/usr/bin/env python
"""Dump best models to disk

Author: Ziga Avsec
Affiliation: TUM
"""
import numpy as np
import pandas as pd
import data
from shutil import copyfile
from glob import glob
from concise.hyopt import CMongoTrials, get_data
import os
from train_all import DIR_ROOT, PROC_DIR, RBP_LIST, RBP_ALL, DB_NAME
from .mongodb_setup import host, port
import logging
import argparse
from pprint import pprint
logging.getLogger().setLevel(logging.WARN)

from joblib import Memory, Parallel, delayed

memory = Memory(cachedir=data.CACHE_DIR, verbose=0)


def print_exp(exp_name, rbp):
    print("-" * 40 + "\nexp: {0}; rbp: {1}".format(exp_name, rbp))


def get_trials(exp_name, rbp):
    return CMongoTrials(DB_NAME, exp_name + "_" + rbp, ip=host, port=port)


def get_models_overall(exp_name, rbp, ignore_cache=False):
    print_exp(exp_name, rbp)

    basedir = "{root}/models/{rbp}/".format(root=DIR_ROOT, rbp=rbp)
    out_h5 = basedir + "{method}.h5".format(method=exp_name)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    trials = get_trials(exp_name, rbp)

    # no trials yet - return None
    if trials.n_ok() == 0:
        return

    # get best trial parameters
    tid = trials.best_trial_tid()
    model_path = trials.get_trial(tid)["result"]["path"]["model"]
    copyfile(model_path, out_h5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run models")
    parser.add_argument("--ignore_cache", action="store_true",
                        help="If enabled, ignore the cached computation")
    parser.add_argument("--n_jobs", type=int, default=10,
                        help="number of jobs to launch in parallel")
    args = parser.parse_args()

    EXPERIMENTS = {"DeepNN_scalar_position_ext_gam": RBP_ALL}

    pprint("Running for experiments:")
    pprint(list(EXPERIMENTS))
    with Parallel(n_jobs=args.n_jobs) as parallel:
        parallel(delayed(get_models_overall)(exp_name, rbp,
                                             ignore_cache=args.ignore_cache)
                 for exp_name, rbps in EXPERIMENTS.items()
                 for rbp in rbps)

"""Post-processing script

- get test predictions
- get loss curves
- get validation accuracies
"""
import os
import matplotlib
import matplotlib.pyplot as plt
from concise.hyopt import CompileFN, CMongoTrials, test_fn
import numpy as np
import pandas as pd
from tqdm import tqdm
from concise.utils.splines import BSpline
from concise.utils.model_data import split_train_test_idx, subset
from keras.models import load_model
import sklearn.metrics as skm
import concise.eval_metrics as cem
import keras.layers as kl
from keras.models import Model
import keras.callbacks as kc
import data
import model
from helper import *
from ..Eclip.precictive_models.mongodb_setup import host, port
import seaborn
import seaborn as sns
import concise.eval_metrics as cem
sns.set(style="whitegrid", color_codes=True)

DB_NAME = "Concise__Splice_branchpoints"
DATA_DIR = os.path.expanduser("~/projects-work/deepcis/data/")
EXP_DIR = DATA_DIR + "/Splice_branchpoints/"

all_trials_names = [["deep_gam", "model_deep"],
                    ["deep_relu", "model_deep_position=relu_scaled"],
                    ["shallow_gam", "model_shallow2"],
                    ["shallow_relu", "model_shallow_position=relu_scaled"],
                    ]

all_trials = [[k, CMongoTrials(DB_NAME, v, ip=host, kill_timeout=30 * 60, port=port)]
              for k, v in all_trials_names]

N_top_tids = 100
for i, (name, t) in enumerate(all_trials):
    print("name: ", name)
    df = t.as_df().sort_values("eval.auprc", ascending=False)
    df.to_csv(EXP_DIR + "/trials/df/{0}.csv".format(name))
    all_trials[i].append({"top_tids": list(df.iloc[:N_top_tids]["tid"])})


N_tids = 1
for i, (name, t, attr) in enumerate(tqdm(all_trials)):
    attr["tid_info"] = {}
    print("name: ", name)
    for tid in attr["top_tids"][:N_tids]:
        print("tid: ", tid)
        param = t.get_param(tid).to_dict()
        if name not in ["shallow_gam_mul", "shallow_gam", "shallow_relu"]:
            # Fix the broken trial names
            param["data"]["truncate"] = False
        train, test = data.data(**param["data"])
        mpath = t.get_trial(tid).to_dict()["result"]["path"]["model"]
        # remove_background_probs(mpath)
        attr["tid_info"][tid] = {"param": param,
                                 "train_all": train,
                                 "data": {"train_valid": train,
                                          "test": test},
                                 "model": load_model(mpath)
                                 }

# save predictions to csv

for name, trials, attr in tqdm(all_trials):
    print("name: ", name)
    for i, tid in enumerate(attr["top_tids"][:N_tids]):
        print("tid: ", tid)
        m = attr["tid_info"][tid]["model"]
        d = attr["tid_info"][tid]["data"]

        y_pred = m.predict(d["test"][0])
        y_true = d["test"][1]
        dt_pred = pd.DataFrame({"y_true": y_true.reshape(
            (-1,)), "y_pred": y_pred.reshape((-1,))})
        dt_pred = dt_pred[dt_pred["y_true"] != -1]
        dt_pred["y_true"] = np.where(dt_pred["y_true"] == 1, "HC", "NEG")
        dt_pred.sort_values("y_pred")
        dt_pred.to_csv(
            EXP_DIR + "/test_predictions/{0}_{1}.csv".format(name, i), index=False)

# get loss history
df_hist = pd.concat([trials.train_history(tid).assign(trial=name)
                     for name, trials, attr in all_trials for tid in attr["top_tids"]])
df_hist.to_csv(EXP_DIR + "/trials/train_history/gam_vs_relu.csv")


# dl = []
# N_tids = 1
# for name, trials, attr in tqdm(all_trials_sub):
#     print("name: ", name)
#     for i, tid in enumerate(attr["top_tids"][:N_tids]):
#         print(EXP_DIR + "/test_predictions/{0}_{1}.csv".format(name, i))
#         df = pd.read_csv(EXP_DIR + "/test_predictions/{0}_{1}.csv".format(name, i))
#         df["name"] = name
#         dl.append(df)

# df = pd.concat(dl)


# df.groupby("name").apply(lambda a: cem.auc(a.y_true != "NEG", a.y_pred))
# df.groupby("name").apply(lambda a: cem.auprc(a.y_true != "NEG", a.y_pred))

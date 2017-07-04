
# coding: utf-8

# In[11]:

import pandas as pd
import numpy as np
from concise.preprocessing import encodeDNA, encodeSplines
from concise.utils.helper import read_json
from concise.hyopt import CMongoTrials, get_data, get_model

import data
import model
from helper import *
DIR_ROOT = "/s/project/deepcis/encode/eclip/"
#DIR_ROOT = "/home/avsec/projects-work/deepcis/data/encode/eclip/"
DB_NAME = "RBP__Eclip"
HOST = "ouga03"

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid", color_codes=True)


# In[2]:

rbp_name = "PUM2"


# In[3]:

# TODO - define 3 different data splits
# 1 existing, 2 new
data_splits = {
    "1_default_split": {"valid_chr": [18, 19, 20, 21, 22, "X"],
                        "test_chr": [1, 2, 3]},
    "2_split": {"valid_chr": [1, 3],
                "test_chr": [2, 4, 6, 8, 10]},
    "3_split": {"valid_chr": [2, 3, 4],
                "test_chr": [5, 6, 7, 8, 9, 10]}
}
split = data_splits["1_default_split"]


# In[4]:

# get existing top hyper-parameters
params = {}
trials = CMongoTrials(DB_NAME, "DeepNN_scalar_position_gam_" + rbp_name, HOST)
params["gam_best"] = trials.get_param(trials.best_trial_tid()).to_dict()
trials = CMongoTrials(DB_NAME, "DeepNN_scalar_position_relu_" + rbp_name, HOST)
params["relu_best"] = trials.get_param(trials.best_trial_tid()).to_dict()


# In[5]:

params


# In[8]:

from keras.callbacks import EarlyStopping


# In[ ]:

dtm_list = []
for k_split in data_splits.keys():
    for k_param in params.keys():
        for ext_pos in ["gam", "relu"]:
            print("{param}, {split}".format(param=k_param, split=k_split))
            param = params[k_param]
            param["model"]["external_pos"]["type"] = ext_pos

            split = data_splits[k_split]
            param["data"] = {**param["data"], **split}
            print("get data...")
            train, valid, test = get_data(data.data, param)
            print("get model...")
            m = get_model(model.model, train, param)
            m.fit(train[0], train[1].astype(int), epochs=150, batch_size=128,
                  validation_data=(valid[0], valid[1].astype(int)),
                  callbacks=[EarlyStopping(patience=5)])
            print("get_metrics...")
            dtm = metrics_dt(m, {"train": (train[0], train[1]),
                                 "valid": (valid[0], valid[1]),
                                 "test": (test[0], test[1])})
            dtm["param"] = k_param
            dtm["split"] = k_split
            dtm["ext_pos"] = ext_pos
            dtm_list.append(dtm)
            print("done!")
dtm_all = pd.concat(dtm_list)
dtm_all.to_csv(DIR_ROOT + "processed/data_split/performances.csv")

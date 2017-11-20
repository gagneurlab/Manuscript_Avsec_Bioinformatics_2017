"""Interpret shallow model
author: Ziga Avsec
"""
import os
from concise.hyopt import CMongoTrials
import numpy as np
import pandas as pd
from concise.utils.model_data import split_train_test_idx, subset
from keras.models import Model
import data
from copy import deepcopy
from helper import plot_pos_dep, logit


def predict_cre(m, trainx):
    inp = m.get_layer("seq").input
    out = m.get_layer("conv_dna_1").output
    nm = Model(inp, out)
    y_pred = nm.predict(trainx)
    m.get_layer("Conv1D_final_layer").get_weights()
    w_seq = m.get_layer("Conv1D_final_layer").get_weights()[0][0][0]
    return y_pred * w_seq


def get_act_kmers(filter_act, filter_len, seqs, lod_thr=0.5):
    # Adopted from Angermueller:
    # https://github.com/cangermueller/deepcpg/blob/6b72aa5a1dc737a7858817e60f6f7365eb5e6470/scripts/dcpg_filter_motifs.py#L114-L151
    assert filter_act.shape[0] == seqs.shape[0]
    idx = np.nonzero(filter_act >= lod_thr)
    kmers = []
    for k in range(len(idx[0])):
        i = int(idx[0][k])
        j = int(idx[1][k])
        kmer = seqs[i, j:(j + filter_len)]
        kmers.append(kmer)
    kmers = np.array(kmers)

    return kmers


if __name__ == "__main__":

    # --------------------------------------------
    # config
    HOST = "ouga03"
    DATA_DIR = os.path.expanduser("~/projects-work/deepcis/data/")
    EXP_NAME = "shallow_spline_trans"
    DB_NAME = "Concise__Splice_branchpoints"
    EXP_DIR = DATA_DIR + "/Splice_branchpoints/"
    # --------------------------------------------
    print("load trials...")
    trial = CMongoTrials(DB_NAME, EXP_NAME, ip=HOST, kill_timeout=30 * 60)
    tid = trial.best_trial_tid()
    m = trial.load_model(tid)
    param = trial.get_param(tid)
    train, test = data.data(**param["data"])
    train_idx, valid_idx = split_train_test_idx(train,
                                                valid_split=.2,
                                                stratified=False, random_state=100)
    valid = subset(train, valid_idx, keep_other=False)
    train_train = subset(train, train_idx)

    y_pred = m.predict(train[0])

    print("extract inferred positional plots...")
    # Extract inferred positional plots
    dfpos = plot_pos_dep(m, param, train, plot=False, use_sigmoid=True)
    # Add some info
    dfpos["type"] = "inferred"

    # Get the measured values
    print("get measured values...")
    dt = pd.read_csv(EXP_DIR + "/processed/branchpointer/train/branchpoint_df_HCN.csv")
    dt = dt.rename(columns={"dist.1": "dist1", "dist.2": "dist2"})
    pos_features = ['dist1', 'dist2',
                    'ppt_start', 'ppt_run_length',
                    'canon_hit1', 'canon_hit2',
                    'canon_hit3', 'canon_hit4', 'canon_hit5']
    # Tidy form
    dtm = pd.melt(dt, id_vars="set", value_vars=pos_features)

    dtm.set = dtm.set == "HC"
    dtmg = dtm.groupby(["variable", "value"]).agg([np.mean, "count"])["set"].reset_index()

    dtmg.rename(columns={"value": "x", "mean": "y", "variable": "feature", "count": "N"}, inplace=True)
    dtmg["type"] = "measured"

    # transform to logit scale
    dtmg.y = logit(dtmg.y) + 2.3

    # Flag and truncate outliers
    dtmg["outlier"] = (dtmg.y < -6) | (dtmg.y > 3)
    dtmg.y[dtmg.y < -6] = -6
    dtmg.y[dtmg.y > 3] = 3

    # Add also the predicted effects
    param2 = deepcopy(param)
    param2["data"]["encode_splines"] = False
    del param2["data"]["n_bases"]
    # Get the data with raw positions
    train, test = data.data(**param2["data"])

    # create long-format df
    df_pos = pd.DataFrame({k: v.reshape((-1)) for k, v in train[0].items() if k in pos_features})
    df_pos["y"] = y_pred.reshape((-1))

    # melt into the appropriate form
    dtm2 = pd.melt(df_pos, id_vars="y", value_vars=pos_features)
    dtmg2 = dtm2.groupby(["variable", "value"]).agg([np.mean, "count"])["y"].reset_index()
    dtmg2.rename(columns={"value": "x", "mean": "y", "variable": "feature", "count": "N"}, inplace=True)
    dtmg2["type"] = "predicted"

    # Back to logit scale
    dtmg2.y = logit(dtmg2.y)
    dtmg2.y = dtmg2.y + 2.3

    # Merge all the features
    dtpos_all = pd.concat([dfpos, dtmg, dtmg2])

    # Restrict ranges
    dtpos_all = dtpos_all[dtpos_all.x < 5000]
    dtpos_all = dtpos_all.groupby("feature", as_index=False).\
        apply(lambda d: d[d.x <= d.x[d.type == "inferred"].max()])
    dtpos_all = dtpos_all.groupby("feature", as_index=False).\
        apply(lambda d: d[d.x >= d.x[d.type == "inferred"].min()])

    # save to csv
    print("save positions to csv...")
    basedir = EXP_DIR + "/interpret/positions/{0}/".format(EXP_NAME)
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    dtpos_all.to_csv(basedir + "/all_positions.csv".format(EXP_NAME))

    # get the PWM's
    # --------------------------------------------
    print("get the pwm's...")
    y_pred_seq = predict_cre(m, train_train[0])

    pwm_model_raw = m.get_layer("conv_dna_1").get_weights()[0][:, :, 0]
    pwm_data = train_train[4][0].pwm

    seqs = get_act_kmers(y_pred_seq, 11, train_train[0]["seq"], lod_thr=2.5)
    pwm_model_transformed = seqs.mean(axis=0)

    # write to text files
    print("save the pwm's...")
    basedir = EXP_DIR + "/interpret/sequence/{0}/".format(EXP_NAME)
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    np.savetxt(basedir + "pwm_model_raw.txt", pwm_model_raw)
    np.savetxt(basedir + "pwm_data.txt", pwm_data)
    np.savetxt(basedir + "pwm_model_transformed.txt", pwm_model_transformed)
    print("Done!")

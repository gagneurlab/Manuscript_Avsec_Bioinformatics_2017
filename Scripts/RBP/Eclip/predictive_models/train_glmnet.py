import data
import os
import pandas as pd
import numpy as np
import gzip
from concise.preprocessing import EncodeSplines
from glmnet import LogitNet
from scipy.sparse import csc_matrix
import concise.eval_metrics as cem
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Imputer
from sklearn.preprocessing import FunctionTransformer
from sklearn.pipeline import make_pipeline

import pickle
from concise.legacy.kmer import kmer_count
from concise.utils.fasta import read_fasta
from concise.utils.helper import write_json
import argparse
from joblib import Parallel, delayed


from train_all import DB_NAME, PROC_DIR, POS_FEATURES, DIR_ROOT, RBP_ALL


from joblib import Memory

memory = Memory(data.CACHE_DIR)


add_eval_metrics = {"auc": cem.auc,
                    "auprc": cem.auprc}

SAVE_ROOT = os.path.join(PROC_DIR, "feature_glmnet_baseline_exp")


def data_kmer(rbp_name,
              n_bases=10,
              kmer=6,
              pos_features=POS_FEATURES,
              valid_chr=[1, 3],
              test_chr=[2, 4, 6, 8, 10]):
    """
    pos_class_weight: positive class weight
    """
    dt_train, dt_valid, dt_test = data.data_split(
        rbp_name + "_extended", valid_chr, test_chr)

    # merge train and valid
    dt_train = pd.concat([dt_train, dt_valid])
    del dt_valid

    seq_train = kmer_count(dt_train.seq.tolist(), kmer)
    seq_test = kmer_count(dt_test.seq.tolist(), kmer)

    # y
    y_train = dt_train.binding_site.as_matrix().reshape(
        (-1, 1)).astype("float")[:, 0]
    y_test = dt_test.binding_site.as_matrix().reshape(
        (-1, 1)).astype("float")[:, 0]

    if n_bases is not None:
        # impute missing values (not part of the pipeline as the Imputer lacks
        # inverse_transform method)
        imp = Imputer(strategy="median")
        imp.fit(dt_train[pos_features])
        dt_train[pos_features] = imp.transform(dt_train[pos_features])
        dt_test[pos_features] = imp.transform(dt_test[pos_features])

        preproc_pipeline = make_pipeline(
            FunctionTransformer(func=data.sign_log_func,
                                inverse_func=data.sign_log_func_inverse))

        # positions
        dtx_train = np.array(dt_train[pos_features])
        dtx_test = np.array(dt_test[pos_features])

        # transform pos features
        preproc_pipeline.fit(dtx_train)
        train_pos = preproc_pipeline.transform(dtx_train)
        test_pos = preproc_pipeline.transform(dtx_test)

        st = EncodeSplines(n_bases=n_bases)

        st.fit(train_pos)

        x_pos_bs_train = st.transform(train_pos)
        x_pos_bs_test = st.transform(test_pos)
        x_pos_bs_train = x_pos_bs_train.reshape((x_pos_bs_train.shape[0], -1))
        x_pos_bs_test = x_pos_bs_test.reshape((x_pos_bs_test.shape[0], -1))

        x_train = np.concatenate([seq_train, x_pos_bs_train], axis=1)
        x_test = np.concatenate([seq_test, x_pos_bs_test], axis=1)
    else:
        st = None
        preproc_pipeline = None
        (x_train, y_train), (x_test, y_test) = (np.array(seq_train),
                                                np.array(y_train)), (np.array(seq_test), np.array(y_test))

    # min-max scale everything
    scaler = MinMaxScaler()
    scaler.fit(x_train)
    x_train = scaler.transform(x_train)
    x_test = scaler.transform(x_test)

    return (x_train, y_train, pos_features, preproc_pipeline), (x_test, y_test)


def train_glmnet(train, test, save_path_pred, save_path_model, save_path_json, n_cores=5):
    ln = LogitNet(alpha=0.5, n_splits=10, n_jobs=n_cores)
    # to sparse
    train_sparse = (csc_matrix(train[0]), csc_matrix(
        train[1].astype(np.float64).reshape((-1, 1))))
    test_sparse = (csc_matrix(test[0]), csc_matrix(
        test[1].astype(np.float64).reshape((-1, 1))))

    print("train the model")
    ln.fit(train_sparse[0], train[1])

    print("get predictions")
    y_pred = ln.predict_proba(test_sparse[0])[:, 1]
    auprc = cem.auprc(test[1], y_pred)
    auc = cem.auc(test[1], y_pred)

    # csv
    print("save csv")
    dt = pd.DataFrame({"y_true": test[1], "y_pred": y_pred})
    dt.to_csv(save_path_pred)

    # json
    print("save json")
    write_json({"auprc": auprc,
                "auc": auc},
               save_path_json)
    # model
    print("save model")
    pickle.dump(ln, open(save_path_model, "wb"))


def train_kmer(protein, n_bases=None, kmer=6, n_cores=5):

    print("Training for protein: {0}".format(protein))
    pos_str = "seq+dist" if n_bases is None else "seq"
    save_path_pred = os.path.join(
        SAVE_ROOT, "predictions", protein, "kmer_k={0}_pos={1}.csv".format(kmer, pos_str))
    save_path_model = os.path.join(
        SAVE_ROOT, "models", protein, "kmer_k={0}_pos={1}.pkl".format(kmer, pos_str))
    save_path_json = os.path.join(
        SAVE_ROOT, "evaluation", protein, "kmer_k={0}_pos={1}.json".format(kmer, pos_str))

    print("save_path_json: " + save_path_json)
    os.makedirs(os.path.dirname(save_path_pred), exist_ok=True)
    os.makedirs(os.path.dirname(save_path_model), exist_ok=True)
    os.makedirs(os.path.dirname(save_path_json), exist_ok=True)

    train, test = data_kmer(protein, n_bases, kmer)

    train_glmnet(train, test, save_path_pred,
                 save_path_model, save_path_json, n_cores)
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run glmnet model")
    parser.add_argument("--n_cores_cv", type=int,
                        default=10, help="Number of cores")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of jobs")
    parser.add_argument("--kmer", type=int, default=6, help="k in k-mer")
    args = parser.parse_args()

    Parallel(n_jobs=args.n_jobs)(delayed(memory.cache(train_kmer))(protein, n_bases, kmer=args.kmer, n_cores=args.n_cores_cv)
                                 for protein in RBP_ALL
                                 for n_bases in [10, None])

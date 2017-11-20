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
import pickle
from concise.legacy.kmer import kmer_count
from concise.utils.fasta import read_fasta
from concise.utils.helper import write_json
import argparse
from joblib import Parallel, delayed
import dill

PROC_DIR = "../../../data/clip/"
from joblib import Memory

memory = Memory(os.path.join(PROC_DIR, "cache"))


def read_seq(protein, split="training"):
    seq_file = os.path.join(data.dtdir, protein, '5000', split + '_sample_0', 'sequences.fa.gz')
    with gzip.open(seq_file, 'r') as fp:
        seq_list = []
        seq = ''
        for line in fp:
            if str(line)[2] == '>':
                if len(seq):
                    seq_list.append(seq.replace("U", "T"))

                seq = ''
            else:
                seq = seq + str(line)[2:-3]
        if len(seq):
            seq_list.append(seq.replace("U", "T"))

    return seq_list

def read_position(protein, split):
    pos_file = os.path.join(data.dtdir, protein, '5000', split + '_sample_0', 'positions.csv')
    return data.read_positions(pos_file)
    
def load_data(protein, n_bases=10, kmer=6):
    path = os.path.join(data.dtdir, protein, '5000')
    
    # sequence
    x_seq_train = kmer_count(read_seq(protein, "training"), kmer)
    x_seq_test = kmer_count(read_seq(protein, "test"), kmer)
  
    # y
    y_train = np.loadtxt(gzip.open(os.path.join(path, 'training_sample_0', "matrix_Response.tab.gz")), skiprows=1)
    y_test = np.loadtxt(gzip.open(os.path.join(path, 'test_sample_0', "matrix_Response.tab.gz")), skiprows=1)
    
    if n_bases is not None:
        # position
        x_pos_train = read_position(protein, "training")
        x_pos_test = read_position(protein, "test")
        
        st = EncodeSplines(n_bases=n_bases)

        st.fit(x_pos_train)

        x_pos_bs_train = st.transform(x_pos_train)
        x_pos_bs_test = st.transform(x_pos_test)

        x_pos_bs_train = x_pos_bs_train.reshape((x_pos_bs_train.shape[0], -1))
        x_pos_bs_test = x_pos_bs_test.reshape((x_pos_bs_test.shape[0], -1))


        x_train = np.concatenate([x_seq_train, x_pos_bs_train], axis = 1)
        x_test = np.concatenate([x_seq_test, x_pos_bs_test], axis = 1)
        #return (x_train, y_train, st), (x_test, y_test)
    else:
        st = None
        (x_train, y_train), (x_test, y_test) = (np.array(x_seq_train), np.array(y_train)), (np.array(x_seq_test), np.array(y_test))
        
    imp = Imputer(missing_values='NaN', strategy='median', axis=0)
    imp.fit(x_train)
    x_train = imp.transform(x_train)
    x_test = imp.transform(x_test)
    
    scaler = MinMaxScaler()
    scaler.fit(x_train)
    x_train = scaler.transform(x_train)
    x_test = scaler.transform(x_test)
    
    return (x_train, y_train, st), (x_test, y_test)


def train_glmnet(train, test, save_path_pred, save_path_model, save_path_json, n_cores=5):
    ln = LogitNet(alpha=0.5, n_splits=10, n_jobs=n_cores)
    # to sparse
    train_sparse = (csc_matrix(train[0]), csc_matrix(train[1].astype(np.float64).reshape((-1,1))))
    test_sparse = (csc_matrix(test[0]), csc_matrix(test[1].astype(np.float64).reshape((-1,1))))
    
    print("train the model")
    ln.fit(train_sparse[0], train[1])
    
    print("get predictions")
    y_pred = ln.predict_proba(test_sparse[0])[:, 1]    
    auprc = cem.auprc(test[1], y_pred)
    auc = cem.auc(test[1], y_pred)
    
    # csv
    print("save csv")
    dt = pd.DataFrame( {"y_true": test[1], "y_pred": y_pred} )
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
    save_path_pred = os.path.join(PROC_DIR, "predictions", protein, "kmer_k={0}_pos={1}.csv".format(kmer, pos_str))
    save_path_model = os.path.join(PROC_DIR, "models", protein, "kmer_k={0}_pos={1}.pkl".format(kmer, pos_str))
    save_path_json = os.path.join(PROC_DIR, "evaluation", protein, "kmer_k={0}_pos={1}.json".format(kmer, pos_str))
    
    os.makedirs(os.path.dirname(save_path_pred), exist_ok=True)
    os.makedirs(os.path.dirname(save_path_model), exist_ok=True)
    os.makedirs(os.path.dirname(save_path_json), exist_ok=True)
    
    train, test = load_data(protein, n_bases, kmer)
    
    train_glmnet(train, test, save_path_pred, save_path_model, save_path_json, n_cores)
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run glmnet model")
    parser.add_argument("--n_cores_cv", type=int, default=10, help="Number of cores")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of jobs")
    parser.add_argument("--kmer", type=int, default=6, help="k in k-mer")
    args = parser.parse_args()
    
    Parallel(n_jobs=args.n_jobs)(delayed(memory.cache(train_kmer))(protein, n_bases, kmer=args.kmer, n_cores=args.n_cores_cv)
                                 for protein in os.listdir(data.dtdir) 
                                 for n_bases in [10, None])

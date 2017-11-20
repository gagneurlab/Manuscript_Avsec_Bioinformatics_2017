import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from concise.preprocessing.sequence import pad_sequences, encodeDNA
from concise.preprocessing.splines import encodeSplines
from concise.preprocessing.structure import encodeRNAStructure
from concise.utils.helper import merge_dicts, read_json, write_json
from concise.utils.pwm import PWM
import os
import pickle
# load the data

# TODO - update
DATA_ROOT = os.path.expanduser("~/projects-work/code_spline_trans/data/")
BR_DATA = DATA_ROOT + "/Splice_branchpoints/processed/branchpointer/"
BR_PWM = DATA_ROOT + "/Splice_branchpoints/processed/branchpointer/hc_pwm.json"
BR_SPLICE_RACK_PATH = DATA_ROOT + "/Splice_branchpoints/raw/U12_splice_sites_SpliceRack/"

pos_columns = ["dist1", "dist2",
               "ppt_start", "ppt_run_length",
               'canon_hit1', 'canon_hit2', 'canon_hit3',
               'canon_hit4', 'canon_hit5']


def truncate_max(x, maxval):
    return np.minimum(x, maxval)


TRUNCATE_VALUES = {
    'dist1': 10000,
    'canon_hit1': 75,
    'canon_hit2': 150,
    'canon_hit3': 150,
    'canon_hit4': 150,
    'canon_hit5': 150,
}


def truncate_values(k, v):
    if k in TRUNCATE_VALUES:
        return truncate_max(v, TRUNCATE_VALUES[k])
    else:
        return v


def get_spliceRack_bp(pseudocount=0.05):
    return [pseudocount + PWM(np.loadtxt(BR_SPLICE_RACK_PATH + "/AT_AC_U12.txt"),
                              "AT_AC_U12"),
            pseudocount + PWM(np.loadtxt(BR_SPLICE_RACK_PATH + "/GT_AG_U12.txt"),
                              "GT_AG_U12")]


def get_branchpoint_pwm_list(cache=True):
    l = []
    if os.path.isfile(BR_PWM) and cache:
        l.append(PWM.from_config(read_json(BR_PWM)))
    else:
        dt = pd.read_csv(DATA_ROOT + "/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")
        # colmeans
        dt = dt[dt.set == "HC"]
        dtseq = dt[dt.columns[dt.columns.str.match("^seq_")]] - 1
        pwm = np.array(dtseq.mean()).reshape((-1, 4))
        assert np.allclose(pwm.sum(1), 1)
        p = PWM(pwm, name="U2_branchpoint")
        write_json(p.get_config(), BR_PWM)
        l.append(p)
    l.append(PWM(0.05 + np.loadtxt(BR_SPLICE_RACK_PATH + "/GT_AG_U12.txt"),
                 "GT_AG_U12_branchpoint"))
    return l


def fill_nan(x, value=0):
    x[np.isnan(x)] = value
    return x


def col2num_array(col, na=np.nan, maxlen=None):
    """Convert a pandas column with comma separated values into numeric np.ndarray
    """
    cons_list = pad_sequences(col.str.split(','), value=["NA"], maxlen=maxlen)
    x_cons = np.asarray(cons_list, dtype='<U20')  # parse strings to have maximal length of 20
    x_cons[x_cons == "NA"] = 'nan'
    x_cons = x_cons.astype(float)
    return x_cons


def get_minmax_range(train_x, test_x):
    train_stats = {k: {"min": min(np.nanmin(train_x[k]), np.nanmin(test_x[k])),
                       "max": max(np.nanmax(train_x[k]), np.nanmax(test_x[k])), }
                   for k in train_x.keys()}
    return train_stats


# NOTE: concise now implements a transoformer API - concise.preprocessing.EncodeSplines with methods:
# .fit
# .predict
# and operates on multiple features in a single array. Please use that one instead.
def encodeSplines_common(x, pos_range, n_bases, spline_order):
    x_processed = {k: fill_nan(encodeSplines(v,
                                             n_bases=n_bases,
                                             spline_order=spline_order,
                                             start=pos_range[k]["min"],
                                             end=pos_range[k]["max"]))
                   for k, v in x.items()}
    return x_processed

def get_structure(seq_vec, cache_path, cache):
    if cache and os.path.isfile(cache_path):
        print("loading cached array")
        return np.load(cache_path)
    else:
        arr = encodeRNAStructure(seq_vec)
        np.save(cache_path, arr)
        return arr


def data(n_bases=10, spline_order=3, pos_class_weight=1.0, truncate=True,
         encode_splines=True, cache=True, minmax_scale=False):
    dtw_train = pd.read_csv(BR_DATA + "/train/wide_data.csv")
    dtw_test = pd.read_csv(BR_DATA + "/test/wide_data.csv")

    # replace some columns names
    dtw_train.rename(columns={"dist.1": "dist1", "dist.2": "dist2"}, inplace=True)
    dtw_test.rename(columns={"dist.1": "dist1", "dist.2": "dist2"}, inplace=True)

    y_train = col2num_array(dtw_train["is_branchpoint"], -1)
    x_train_pos = {col: col2num_array(dtw_train[col]) for col in pos_columns}
    x_train_seq = encodeDNA(dtw_train["seq"])
    y_test = col2num_array(dtw_test["is_branchpoint"])
    x_test_pos = {col: col2num_array(dtw_test[col]) for col in pos_columns}
    x_test_seq = encodeDNA(dtw_test["seq"])

    # get secondary structure - use a larger sequence window to precompute it
    # to get the larger sequence window first
    assert dtw_train["seq_wide"][0][100:-100] == dtw_train["seq"][0]  # correct values
    x_train_struct = get_structure(dtw_train["seq_wide"], BR_DATA + "/train/rna_structure.npy", cache)
    x_train_struct = x_train_struct[:, 100:-100]
    x_test_struct = get_structure(dtw_test["seq_wide"], BR_DATA + "/test/rna_structure.npy", cache)
    x_test_struct = x_train_struct[:, 100:-100]

    # convert all np.nans into -1's for y's and 0's for x's
    if truncate:
        x_train_pos = {k: truncate_values(k, v) for k, v in x_train_pos.items()}
        x_test_pos = {k: truncate_values(k, v) for k, v in x_test_pos.items()}

    # Ranges computed only on the trainin set
    position_stats = get_minmax_range(x_train_pos, x_train_pos)

    minmax_scalers = None
    if encode_splines:
        x_train_pos_spl = encodeSplines_common(x_train_pos, position_stats, n_bases, spline_order)
        x_test_pos_spl = encodeSplines_common(x_test_pos, position_stats, n_bases, spline_order)
    else:
        x_train_pos_spl = {k: fill_nan(v)[:, :, np.newaxis] for k, v in x_train_pos.items()}
        x_test_pos_spl = {k: fill_nan(v)[:, :, np.newaxis] for k, v in x_test_pos.items()}
        if minmax_scale:
            minmax_scalers = {k: MinMaxScaler().fit(v.reshape(-1, 1)) for k, v in x_train_pos_spl.items()}
            x_train_pos_spl = {k: minmax_scalers[k].transform(v.reshape((-1, 1))).reshape(v.shape)
                               for k, v in x_train_pos_spl.items()}
            x_test_pos_spl = {k: minmax_scalers[k].transform(v.reshape((-1, 1))).reshape(v.shape)
                              for k, v in x_test_pos_spl.items()}

    x_train = merge_dicts(x_train_pos_spl, {"seq": x_train_seq, "struct": x_train_struct})
    x_test = merge_dicts(x_test_pos_spl, {"seq": x_test_seq, "struct": x_test_struct})

    y_train = fill_nan(y_train, -1)[:, :, np.newaxis]
    y_test = fill_nan(y_test, -1)[:, :, np.newaxis]

    # todo update return
    sample_weight = np.squeeze(np.where(y_train == 1, pos_class_weight, 1), -1)
    return (x_train, y_train, sample_weight,
            pos_columns, get_branchpoint_pwm_list(), position_stats, minmax_scalers), (x_test, y_test)

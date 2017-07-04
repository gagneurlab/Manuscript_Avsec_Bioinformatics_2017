import pandas as pd
import numpy as np
from concise.preprocessing import encodeDNA, encodeSplines
from concise.utils.helper import read_json
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import FunctionTransformer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Imputer
from tempfile import mkdtemp
from joblib import Memory

DIR_ROOT = "/s/project/deepcis/encode/eclip/"
# DIR_ROOT = "/home/avsec/projects-work/deepcis/data/encode/eclip/"
CACHE_DIR = DIR_ROOT + "cache/"

memory = Memory(cachedir=CACHE_DIR, verbose=0)


def get_pos_ranges(d):
    min_val = min([v.min() for k, v in d.items()])
    max_val = max([v.max() for k, v in d.items()])
    return {"min": min_val, "max": max_val}


rbp_name = "UPF1"


def data_split(rbp_name, valid_chr=[1, 3],
               test_chr=[2, 4, 6, 8, 10]):

    assert len(set(valid_chr).intersection(test_chr)) == 0
    valid_chr = ["chr" + str(i) for i in valid_chr]
    test_chr = ["chr" + str(i) for i in test_chr]

    dt_train = pd.read_csv(DIR_ROOT + "/processed/design_matrix/train/{0}.csv".
                           format(rbp_name))
    dt_valid = pd.read_csv(DIR_ROOT + "/processed/design_matrix/valid/{0}.csv".
                           format(rbp_name))
    dt_test = pd.read_csv(DIR_ROOT + "/processed/design_matrix/test/{0}.csv".
                          format(rbp_name))

    dt_train.rename(columns={"TSS": "tss"}, inplace=True)
    dt_valid.rename(columns={"TSS": "tss"}, inplace=True)
    dt_test.rename(columns={"TSS": "tss"}, inplace=True)
    dt_all = pd.concat([dt_train, dt_valid, dt_test])
    dt_valid = dt_all[dt_all.seqnames.isin(valid_chr)]
    dt_test = dt_all[dt_all.seqnames.isin(test_chr)]
    dt_train = dt_all[~dt_all.seqnames.isin(valid_chr + test_chr)]
    return dt_train, dt_valid, dt_test

def calc_perc(data):
    train, valid, test = data

    n_train = len(train)
    n_valid = len(valid)
    n_test = len(test)
    return np.array([n_train, n_valid, n_test]) / (n_train + n_valid + n_test)


POS_FEATURES = ['tss', 'polya', 'exon_intron', 'intron_exon', 'start_codon',
                'stop_codon', 'gene_start', 'gene_end']


def sign_log_func(x):
    return np.sign(x) * np.log10(np.abs(x) + 1)


def sign_log_func_inverse(x):
    return np.sign(x) * (np.power(10, np.abs(x)) - 1)


# memoize the function
# @memory.cache
def data_extended(rbp_name, n_bases=10,
                  pos_class_weight=1.0,
                  scale="sign_log",  # or "nat"
                  valid_chr=[1, 3],
                  test_chr=[2, 4, 6, 8, 10]):
    """
    pos_class_weight: positive class weight
    """
    dt_train, dt_valid, dt_test = data_split(rbp_name + "_extended", valid_chr, test_chr)

    seq_train = encodeDNA(dt_train.seq.tolist())
    seq_valid = encodeDNA(dt_valid.seq.tolist())
    seq_test = encodeDNA(dt_test.seq.tolist())

    # impute missing values (not part of the pipeline as the Imputer lacks inverse_transform method)
    imp = Imputer(strategy="median")
    imp.fit(pd.concat([dt_train[POS_FEATURES], dt_valid[POS_FEATURES]]))
    dt_train[POS_FEATURES] = imp.transform(dt_train[POS_FEATURES])
    dt_valid[POS_FEATURES] = imp.transform(dt_valid[POS_FEATURES])
    dt_test[POS_FEATURES] = imp.transform(dt_test[POS_FEATURES])

    if scale == "sign_log":
        preproc_pipeline = make_pipeline(
            FunctionTransformer(func=sign_log_func,
                                inverse_func=sign_log_func_inverse),
            MinMaxScaler()
        )
    elif scale == "nat":
        preproc_pipeline = make_pipeline(
            MinMaxScaler()
        )
    else:
        ValueError("scale argument invalid")

    # transform pos features
    preproc_pipeline.fit(pd.concat([dt_train[POS_FEATURES], dt_valid[POS_FEATURES]]))
    train_pos = pd.DataFrame(preproc_pipeline.transform(dt_train[POS_FEATURES]), columns=POS_FEATURES)
    valid_pos = pd.DataFrame(preproc_pipeline.transform(dt_valid[POS_FEATURES]), columns=POS_FEATURES)
    test_pos = pd.DataFrame(preproc_pipeline.transform(dt_test[POS_FEATURES]), columns=POS_FEATURES)

    def create_feature_dict(dt):
        raw_dist = {"raw_dist_" + k: np.array(v)[:, np.newaxis, np.newaxis] for k, v in dt.items()}
        dist = {"dist_" + k: encodeSplines(np.array(v)[:, np.newaxis], start=0, end=1) for k, v in dt.items()}
        return {**raw_dist, **dist}

    train_dist = create_feature_dict(train_pos)
    valid_dist = create_feature_dict(valid_pos)
    test_dist = create_feature_dict(test_pos)

    x_train = {"seq": seq_train, **train_dist}
    x_valid = {"seq": seq_valid, **valid_dist}
    x_test = {"seq": seq_test, **test_dist}

    # y
    y_train = dt_train.binding_site.as_matrix().reshape((-1, 1)).astype("float")
    y_valid = dt_valid.binding_site.as_matrix().reshape((-1, 1)).astype("float")
    y_test = dt_test.binding_site.as_matrix().reshape((-1, 1)).astype("float")
    sample_weight = np.squeeze(np.where(y_train == 1, pos_class_weight, 1), -1)

    return (x_train, y_train, sample_weight, POS_FEATURES, preproc_pipeline), \
        (x_valid, y_valid),\
        (x_test, y_test)


# @memory.cache
def data(rbp_name, n_bases=10,
         pos_class_weight=1.0,
         tss_trunc=2000, polya_trunc=2000,
         pos_as_track=False, kernel_size=10,
         scale_raw=False,
         valid_chr=[1, 3],
         test_chr=[2, 4, 6, 8, 10]):
    """
    pos_class_weight: positive class weight
    """
    # info = read_json(DIR_ROOT + "/processed/design_matrix/meta_info/{0}.json".
    #                  format(rbp_name))

    dt_train, dt_valid, dt_test = data_split(rbp_name, valid_chr, test_chr)

    # TODO - not working just with dt_train.seq ?!?!?
    seq_train = encodeDNA(dt_train.seq.tolist())
    seq_valid = encodeDNA(dt_valid.seq.tolist())
    seq_test = encodeDNA(dt_test.seq.tolist())

    tss_dist = {"train": dt_train.TSS_distance.values,
                "valid": dt_valid.TSS_distance.values,
                "test": dt_test.TSS_distance.values}
    polya_dist = {"train": dt_train.polya_distance.values,
                  "valid": dt_valid.polya_distance.values,
                  "test": dt_test.polya_distance.values}

    seq_length = seq_train.shape[1]
    pos_length = seq_length - kernel_size + 1

    def expand_positions(x, pos_length):
        """If pos_as_track, use it"""
        x = x.reshape((-1, 1))
        # 1. create a matrix with incrementing positions
        incr_array = np.arange(pos_length) - pos_length // 2
        # expand to have the same shape as x
        positions_offset = np.repeat(incr_array.reshape((1, -1)), x.shape[0], axis=0)
        return positions_offset + x

    if pos_as_track:
        tss_dist = {k: expand_positions(v, pos_length) for k, v in tss_dist.items()}
        polya_dist = {k: expand_positions(v, pos_length) for k, v in polya_dist.items()}
        shift = pos_length // 2 + 2
    else:
        tss_dist = {k: v[:, np.newaxis] for k, v in tss_dist.items()}
        polya_dist = {k: v[:, np.newaxis] for k, v in polya_dist.items()}
        shift = 1

    # transform polya_distance - change order
    tss_dist = {k: (v + shift) for k, v in tss_dist.items()}
    polya_dist = {k: -1 * (v - shift) for k, v in polya_dist.items()}

    tss_pos_ranges = get_pos_ranges(tss_dist)
    polya_pos_ranges = get_pos_ranges(polya_dist)

    def get_tss_nat_dist(x):
        return encodeSplines(x, n_bases=n_bases,
                             start=tss_pos_ranges["min"],
                             end=tss_trunc)

    def get_tss_log_dist(x):
        return encodeSplines(np.log10(x), n_bases=n_bases,
                             start=np.log10(tss_pos_ranges["min"]),
                             end=np.log10(tss_pos_ranges["max"]),
                             )

    def get_polya_nat_dist(x):
        return encodeSplines(x, n_bases=n_bases,
                             start=polya_pos_ranges["min"],
                             end=polya_trunc)

    def get_polya_log_dist(x):
        return encodeSplines(np.log10(x), n_bases=n_bases,
                             start=np.log10(polya_pos_ranges["min"]),
                             end=np.log10(polya_pos_ranges["max"]),
                             )

    # min-max scaler
    mms_tss = MinMaxScaler()
    mms_tss.fit(np.log10(tss_dist["train"]).reshape((-1, 1)))
    mms_polya = MinMaxScaler()
    mms_polya.fit(np.log10(polya_dist["train"]).reshape((-1, 1)))

    def get_raw_tss_log_dist(x):
        sh = x.shape
        if scale_raw:
            return mms_tss.transform(np.log10(x).reshape((-1, 1))).\
                reshape(sh)[:, :, np.newaxis]
        else:
            return np.log10(x)[:, :, np.newaxis]

    def get_raw_polya_log_dist(x):
        sh = x.shape
        if scale_raw:
            return mms_polya.transform(np.log10(x).reshape((-1, 1))).\
                reshape(sh)[:, :, np.newaxis]
        else:
            return np.log10(x)[:, :, np.newaxis]

    y_train = dt_train.binding_site.as_matrix().reshape((-1, 1)).astype("float")
    y_valid = dt_valid.binding_site.as_matrix().reshape((-1, 1)).astype("float")
    y_test = dt_test.binding_site.as_matrix().reshape((-1, 1)).astype("float")
    sample_weight = np.squeeze(np.where(y_train == 1, pos_class_weight, 1), -1)
    return ({"seq": seq_train,
             "dist_tss_nat": get_tss_nat_dist(tss_dist["train"]),
             "dist_tss_log": get_tss_log_dist(tss_dist["train"]),
             "dist_polya_nat": get_polya_nat_dist(polya_dist["train"]),
             "dist_polya_log": get_polya_log_dist(polya_dist["train"]),
             # "raw_dist_tss_nat": tss_dist["train"], # Not supported, not thresholding it
             "raw_dist_tss_log": get_raw_tss_log_dist(tss_dist["train"]),
             # "raw_dist_polya_nat": polya_dist["train"],
             "raw_dist_polya_log": get_raw_polya_log_dist(polya_dist["train"])},
            y_train, sample_weight, tss_pos_ranges, polya_pos_ranges, mms_tss, mms_polya),\
        ({"seq": seq_valid,
          "dist_tss_nat": get_tss_nat_dist(tss_dist["valid"]),
          "dist_tss_log": get_tss_log_dist(tss_dist["valid"]),
          "dist_polya_nat": get_polya_nat_dist(polya_dist["valid"]),
          "dist_polya_log": get_polya_log_dist(polya_dist["valid"]),
          # "raw_dist_tss_nat": tss_dist["valid"],
          "raw_dist_tss_log": get_raw_tss_log_dist(tss_dist["valid"]),
          # "raw_dist_polya_nat": polya_dist["valid"],
          "raw_dist_polya_log": get_raw_polya_log_dist(polya_dist["valid"])},
         y_valid),\
        ({"seq": seq_test,
          "dist_tss_nat": get_tss_nat_dist(tss_dist["test"]),
          "dist_tss_log": get_tss_log_dist(tss_dist["test"]),
          "dist_polya_nat": get_polya_nat_dist(polya_dist["test"]),
          "dist_polya_log": get_polya_log_dist(polya_dist["test"]),
          # "raw_dist_tss_nat": tss_dist["test"],
          "raw_dist_tss_log": get_raw_tss_log_dist(tss_dist["test"]),
          # "raw_dist_polya_nat": polya_dist["test"],
          "raw_dist_polya_log": get_raw_polya_log_dist(polya_dist["test"])},
         y_test)

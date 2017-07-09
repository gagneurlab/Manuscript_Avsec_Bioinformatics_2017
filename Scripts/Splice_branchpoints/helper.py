import sklearn.metrics as skm
import concise.eval_metrics as cem
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pprint import pprint
from concise.utils.splines import BSpline
from concise.hyopt import eval_model
import concise.eval_metrics as cem
import keras.layers as kl
from keras.models import Model

def plot_pr_curve(y_true, y_pred, show=True):
    y_true, y_pred = cem._mask_value(y_true, y_pred)
    auprc_macro = skm.average_precision_score(y_true, y_pred, average="macro")
    #auprc_micro = skm.average_precision_score(y_true, y_pred, average="micro")
    precision, recall, _ = skm.precision_recall_curve(y_true, y_pred)
    plt.plot(recall, precision, label='Precision-Recall curve')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall: AUC={0:0.2f}'.format(auprc_macro))
    if show:
        plt.show()

def plot_roc_curve(y_true, y_pred, show=True):
    y_true, y_pred = cem._mask_value(y_true, y_pred)
    fpr, tpr, _ = skm.roc_curve(y_true, y_pred)
    roc_auc = skm.auc(fpr, tpr)
    plt.plot(fpr, tpr)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('ROC: AUC={0:0.2f}'.format(roc_auc))
    if show:
        plt.show()


def param2str(param):
    from concise.hyopt import _dict_to_filestring, _flatten_dict_ignore
    """represent model with a short string"""
    keep = ["n_bases", "filters", "nonlinearity", "truncate", "pos_cls_w",
            "use_pssm", "type", "n_hidden", "lr", "dropout_rate", "n_layers"]
    import uuid

    rnd = str(uuid.uuid4())[:4]
    flat_dict = _flatten_dict_ignore(param)
    fs_str = _dict_to_filestring({k: flat_dict.get(k, None) for k in keep})
    fs_str = fs_str.replace("nonlinearity", "nl")
    fs_str = fs_str.replace("n_hidden", "hid_n")
    fs_str = fs_str.replace("n_layers", "hid_nl")
    fs_str = fs_str.replace("filters", "filt")
    fs_str = fs_str.replace("concatenate", "concat")
    fs_str = fs_str.replace("dropout_rate", "hid_d")
    fs_str = fs_str.replace("truncate", "trunc")
    return fs_str + ";" + rnd

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def logit(x):
    return np.log(x / (1 - x))


# TODO - put this plot into trials object
def plot_history(trials, tid, scores=["loss", "f1", "accuracy"],
                 figsize=(15, 3)):
    """Plot the loss curves"""
    history = trials.train_history(tid)
    fig = plt.figure(figsize=figsize)
    for i, score in enumerate(scores):
        plt.subplot(1, len(scores), i + 1)

        plt.plot(history[score], label="train")
        plt.plot(history['val_' + score], label="validation")
        plt.tight_layout()
        plt.title(score)
        plt.ylabel(score)
        plt.xlabel('epoch')
        plt.legend(loc='best')
    return fig


def plot_pos_dep(m, param, train, use_sigmoid=False, plot=True, figsize=(15, 6)):
    position_stats = train[5]
    # pprint(position_stats)
    w_final = m.layers[-1].get_weights()[0][0, :]
    #print("w_final: ", w_final)
    df_list = []
    if plot:
        plt.figure(figsize=figsize)
    for i, l in enumerate(train[3]):
        # Get the range
        start = position_stats[l]["min"]
        end = position_stats[l]["max"]
        bs = BSpline(start, end, n_bases=param["data"]["n_bases"])
        x_range = np.linspace(start, end, 1000)
        X_pred = bs.predict(x_range)

        if param["model"]["pos_effect"]["merge"]["type"] == "concatenate":
            w_cur = w_final.reshape((-1, ))[1 + i]
        else:
            w_cur = w_final.reshape((1, -1))

        layer = m.get_layer(l + "_conv")
        w = layer.get_weights()[0].squeeze(0)
        y_pred = np.dot(X_pred, w) * w_cur
        if plot:
            plt.subplot(3, 3, i + 1)
            plt.tight_layout()
            if use_sigmoid:
                plt.plot(x_range, sigmoid(y_pred))
            else:
                plt.plot(x_range, y_pred)
            plt.xlim(start, end)
            plt.title(l)
        y = y_pred.reshape((-1,))
        assert len(x_range) == len(y)
        df_list.append(pd.DataFrame({"x": x_range, "y": y, "feature": l}))

    dfpos = pd.concat(df_list)
    return dfpos

def plot_pos_dep_relu(m, param, train, scaled=False, use_sigmoid=False, figsize=(15, 6)):
    position_stats = train[5]
    # pprint(position_stats)
    w_final = m.layers[-1].get_weights()[0][0, :]
    df_list = []
    plt.figure(figsize=figsize)
    for i, l in enumerate(train[3]):
        # Get the range
        input_l = kl.Input((1000, 1), name=l)
        out1 = m.get_layer(l + "_conv_features")(input_l)
        out = m.get_layer(l + "_conv_combine")(out1)
        cm = Model(input_l, out)

        start = position_stats[l]["min"]
        end = position_stats[l]["max"]
        if scaled:
            x_range = np.linspace(0, 1, 1000)
        else:
            x_range = np.linspace(start, end, 1000)

        x_range = x_range.reshape((1, -1, 1))
        y = cm.predict(x_range).reshape((-1))

        if scaled:
            x_range = train[6][l].inverse_transform(x_range.reshape((-1, 1)))

        if param["model"]["pos_effect"]["merge"]["type"] == "concatenate":
            w_cur = w_final.reshape((-1, ))[1 + i]
        else:
            w_cur = w_final.reshape((1, -1))

        y_pred = y * w_cur
        plt.subplot(3, 3, i + 1)
        plt.tight_layout()
        x_range = x_range.reshape((-1))
        if use_sigmoid:
            plt.plot(x_range, sigmoid(y_pred))
        else:
            plt.plot(x_range, y_pred)
        plt.xlim(start, end)
        plt.title(l)
        y = y_pred.reshape((-1,))
        assert len(x_range) == len(y)
        df_list.append(pd.DataFrame({"x": x_range, "y": y, "feature": l}))

    dfpos = pd.concat(df_list)
    return dfpos


def plot_roc_pr(m, dataset_list, labels=[], figsize=(6, 2.5), return_data=False):
    data = [(d[1], m.predict(d[0])) for d in dataset_list]
    plt.figure(figsize=figsize)
    plt.subplot(121)
    for y_true, y_pred in data:
        plot_roc_curve(y_true, y_pred, show=False)
    plt.legend(labels, loc='best')
    plt.tight_layout()
    plt.subplot(122)
    for y_true, y_pred in data:
        plot_pr_curve(y_true, y_pred, show=False)
    plt.legend(labels, loc='best')
    plt.tight_layout()
    if return_data:
        return data


def metrics_dt(m, datasets, add_eval_metrics={"auc": cem.auc, "auprc": cem.auprc}):
    """Return a table with model metrics as columns"""
    data = [{"dataset": k, **eval_model(m, d, add_eval_metrics)} for k, d in datasets.items()]
    colorder = ["dataset"] + m.metrics_names + list(add_eval_metrics.keys())
    return pd.DataFrame(data)[colorder]

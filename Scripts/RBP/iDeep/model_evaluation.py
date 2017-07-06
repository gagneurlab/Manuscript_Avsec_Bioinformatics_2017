"""Compute the performance of the 3 models for each experiment

Author: Mohammadamin Barekatain
Affiliation: TUM
"""

import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
np.random.seed(1337)

dtdir = '../../../Data/RBP/iDeep/clip/'
outdir = "../../../Output/RBP/iDeep/"
ideep_original = os.path.join(outdir, 'iDeep_original')
ideep_positions_nat = os.path.join(outdir, 'iDeep_scaler_position_relu')
ideep_positions_gam = os.path.join(outdir, 'iDeep_scaler_position_gam')
results_type = np.dtype([('protein', np.str_, 128), ('y_true', np.int0, (1000,)), ('y_pred_proba', np.float64, (1000,)), ('auc', np.float64)])


# I changed 'Mut FUS' to 'Mut-FUS' in the proteinnames.txt so that the delimiter can be ' '
def read_protein_name(filename='proteinnames.txt'):
    protein_dict = {}
    with open(filename, 'r') as fp:
        for line in fp:
            values = line.split(' ')
            key_name = values[0]
            if len(values) == 2:
                protein_dict[key_name] = values[1][:-1]
            else:
                protein_dict[key_name] = values[1] + ' ' + values[2][:-1]
    return protein_dict


def load_result(path):
    protein_dict = read_protein_name()
    results_list = []
    for protein in os.listdir(dtdir):
        if protein == '.DS_Store':
            continue
        dt = pd.read_csv(os.path.join(path, protein + ".csv"))
        auc = roc_auc_score(dt["y_true"], dt["y_pred_proba"])
        protein = protein_dict[protein]
        results_list.append((protein, dt["y_true"], dt["y_pred_proba"], auc))
    return np.array(results_list, dtype=results_type)


def make_all_dt(results_list, method_list, measure, num_samples=200):

    def bootstrap(y_true, y_pred, num_samples, measure):
        n = len(y_true)
        idx = np.random.randint(0, n, (num_samples, n))
        return np.array([measure(y_true[idx[i]], y_pred[idx[i]])
                         for i in range(idx.shape[0])])

    dt_list = []
    for result, method in zip(results_list, method_list):
        for row in result:
            dt = pd.DataFrame()
            dt['rbp'] = [row['protein']] * num_samples
            dt['method'] = [method] * num_samples
            if measure == roc_auc_score:
                dt['auc'] = bootstrap(row['y_true'], row['y_pred_proba'], num_samples, measure)
            elif measure == average_precision_score:
                dt['auprc'] = bootstrap(row['y_true'], row['y_pred_proba'], num_samples, measure)
            dt_list.append(dt)

    return pd.concat(dt_list)


ideep = load_result(ideep_original)
ideep = np.sort(ideep, order='auc')

ideep_positions_nat = load_result(ideep_positions_nat)
ideep_positions_nat = np.sort(ideep_positions_nat, order='auc')

ideep_positions_gam = load_result(ideep_positions_gam)
ideep_positions_gam = np.sort(ideep_positions_gam, order='auc')


dtb_auprc = make_all_dt([ideep, ideep_positions_gam, ideep_positions_nat],
                        ['iDeep', 'iDeep_scaler_position_gam', 'iDeep_scaler_position_relu'], average_precision_score)
dtb_auprc.to_csv(os.path.join(outdir, "iDeep_auprc.csv"))

dtb_auc = make_all_dt([ideep, ideep_positions_gam, ideep_positions_nat],
                      ['iDeep', 'iDeep_scaler_position_gam', 'iDeep_scaler_position_relu'], roc_auc_score)
dtb_auc.to_csv(os.path.join(outdir, "iDeep_auc.csv"))

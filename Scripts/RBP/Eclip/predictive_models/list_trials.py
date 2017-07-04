from concise.hyopt import CMongoTrials
import numpy as np
from glob import glob
import os

exp_name = "DeepNN_scalar_position_gam_all"
# PROC_DIR = "/s/project/deepcis/encode/eclip/processed"
PROC_DIR = "/home/avsec/projects-work/deepcis/data/encode/eclip/processed"
RBP_ALL = [os.path.basename(x).replace(".csv", "")
           for x in glob(PROC_DIR + "/design_matrix/train/*.csv")]

DB_NAME = "RBP__Eclip"

lengths = {}
for rbp in RBP_ALL:
    # config
    c_exp_name = exp_name + "_" + rbp

    # run
    trials = CMongoTrials(DB_NAME, c_exp_name, ip='localhost', kill_timeout=30 * 60)

    print("exp: {exp}, len: {l}".format(exp=c_exp_name, l=len(trials)))
    lengths[c_exp_name] = len(trials)

    # Almost all went through

# element counts
unique, counts = np.unique(np.array([x for x in lengths.values()]),
                           return_counts=True)
print(np.asarray((unique, counts)).T)

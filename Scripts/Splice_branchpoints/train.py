#!/usr/bin/env python
"""Train the models
"""
from hyperopt import fmin, tpe, hp, pyll
from concise.hyopt import CompileFN, CMongoTrials, test_fn
from copy import deepcopy
import numpy as np
from ..Eclip.precictive_models.mongodb_setup import host, port
import data
import model


def print_exp(exp_name):
    print("-" * 40 + "\nexp_name: " + exp_name)


KILL_TIMEOUT = 60 * 20  # 20 minutes

DB_NAME = "Concise__Splice_branchpoints"
# --------------------------------------------
exp_name = "model_deep"
print_exp(exp_name)
# -----
fn = CompileFN(DB_NAME, exp_name,
               data_fn=data.data,
               model_fn=model.model,
               add_eval_metrics=["auprc", "auc"],
               loss_metric="auprc",
               loss_metric_mode="max",
               use_tensorboard=False,
               valid_split=0.2,
               random_state=100)

hyper_params = {
    "data": {"n_bases": 10 + 10 * hp.randint("n_bases", 3),
             "pos_class_weight": 2.0,
             "truncate": hp.choice("truncate", (False, True)),
             # **NOTE**: there was a bug at the time I ran tris experiment. We were actually using truncate=False all the time.
             },
    "model": {"nonlinearity": "relu",
              "filters": 15 * 2**hp.randint("filters", 5),
              "init_motifs": None,
              "pos_effect": {"l2_smooth": hp.loguniform("m_pos_effect_l2_smooth", np.log(1e-18), np.log(1e-1)),
                             "l2": 0,
                             "use_bias": False,
                             "merge": {
                                 "type": "concatenate",
                                 "hidden_fc": {"n_hidden": 10 + 30 * hp.randint("n_hidden", 3),
                                               "activation": "relu",
                                               "dropout_rate": hp.uniform("m_dropout_rate", 0, 0.6),
                                               "n_layers": hp.randint("m_n_layers", 4),
                                               },
              },
              },
              "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(0.005)),
              },
    "fit": {"epochs": 150,
            "patience": 4,
            "batch_size": 128,
            }
}

# test_fn(fn, hyper_params, n_train=1000)
trials = CMongoTrials(DB_NAME, ip=host, exp_name, kill_timeout=30 * 60)
best = fmin(fn, hyper_params, trials=trials, algo=tpe.suggest, max_evals=1000)
print("best_parameters: " + str(best))

# ----------------
exp_name = "model_shallow2"
print_exp(exp_name)
# -----
fn = CompileFN(DB_NAME, exp_name,
               data_fn=data.data,
               model_fn=model.model,
               add_eval_metrics=["auprc", "auc"],
               loss_metric="auprc",
               loss_metric_mode="max",
               use_tensorboard=False,
               valid_split=0.2,
               random_state=100)

hyper_params = {
    "data": {"n_bases": 10 + 10 * hp.randint("n_bases", 3),
             "pos_class_weight": 2.0,
             "truncate": True,
             },
    "model": {"nonlinearity": "relu",
              "filters": 1,
              "init_motifs": {"use_pssm": hp.choice("m_use_pssm", (False, True)),
                              "stddev": hp.loguniform("m_stddev", np.log(1e-5), np.log(0.3)),
                              },
              "pos_effect": {"l2_smooth": hp.loguniform("m_pos_effect_l2_smooth",
                                                        np.log(1e-10),
                                                        np.log(1e2)),
                             "l2": hp.loguniform("m_pos_effect_l2",
                                                 np.log(1e-18),
                                                 np.log(1e2)),
                             "use_bias": False,
                             "merge": hp.choice("m_merge", (
                                 {"type": "multiply"},
                                 {"type": "add"},
                                 {"type": "concatenate", "hidden_fc": None})),
                             },
              "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(0.005)),
              },
    "fit": {"epochs": 150,
            "patience": 2,
            "batch_size": 128,
            }
}

# test_fn(fn, hyper_params, n_train=1000)
trials = CMongoTrials(DB_NAME, ip=host, exp_name, kill_timeout=30 * 60)
best = fmin(fn, hyper_params, trials=trials, algo=tpe.suggest, max_evals=200)
print("best_parameters: " + str(best))

# ----------------
exp_name = "model_shallow_position=relu"
print_exp(exp_name)
# -----
fn = CompileFN(DB_NAME, exp_name,
               data_fn=data.data,
               model_fn=model.model,
               add_eval_metrics=["auprc", "auc"],
               loss_metric="auprc",
               loss_metric_mode="max",
               use_tensorboard=False,
               valid_split=0.2,
               random_state=100)

hyper_params = {
    "data": {"pos_class_weight": 2.0,
             "truncate": hp.choice("d_truncate", (False, True)),
             "encode_splines": False,
             },
    "model": {"nonlinearity": "relu",
              "filters": 1,
              "init_motifs": {"use_pssm": hp.choice("m_use_pssm", (False, True)),
                              "stddev": hp.loguniform("m_stddev", np.log(1e-5), np.log(0.3)),
                              },
              "pos_effect": {"type": "relu",
                             "n_bases": 10 + 10 * hp.randint("n_bases", 3),
                             "l2": hp.loguniform("m_pos_effect_l2",
                                                 np.log(1e-18),
                                                 np.log(1e2)),
                             "use_bias": hp.choice("m_pos_use_bias", (False, True)),
                             "activation": hp.choice("m_pos_activation", (None, "relu")),
                             "merge": hp.choice("m_merge", (
                                 {"type": "multiply"},
                                 {"type": "add"},
                                 {"type": "concatenate", "hidden_fc": None})),
                             },
              "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(0.005)),
              },
    "fit": {"epochs": 150,
            "patience": 2,
            "batch_size": 128,
            }
}

# test_fn(fn, hyper_params, n_train=1000)
trials = CMongoTrials(DB_NAME, ip=host, exp_name, kill_timeout=30 * 60)
best = fmin(fn, hyper_params, trials=trials, algo=tpe.suggest, max_evals=200)
print("best_parameters: " + str(best))
# ----------------
exp_name = "model_shallow_position=relu_scaled"
print_exp(exp_name)
# -----
fn = CompileFN(DB_NAME, exp_name,
               data_fn=data.data,
               model_fn=model.model,
               add_eval_metrics=["auprc", "auc"],
               loss_metric="auprc",
               loss_metric_mode="max",
               use_tensorboard=False,
               valid_split=0.2,
               random_state=100)

hyper_params = {
    "data": {"pos_class_weight": 2.0,
             "truncate": hp.choice("d_truncate", (False, True)),
             "encode_splines": False,
             "minmax_scale": True,
             },
    "model": {"nonlinearity": "relu",
              "filters": 1,
              "init_motifs": {"use_pssm": hp.choice("m_use_pssm", (False, True)),
                              "stddev": hp.loguniform("m_stddev", np.log(1e-5), np.log(0.3)),
                              },
              "pos_effect": {"type": "relu",
                             "n_bases": 10 + 10 * hp.randint("n_bases", 3),
                             "l2": hp.loguniform("m_pos_effect_l2",
                                                 np.log(1e-18),
                                                 np.log(1e2)),
                             "use_bias": hp.choice("m_pos_use_bias", (False, True)),
                             "activation": hp.choice("m_pos_activation", (None, "relu")),
                             "merge": hp.choice("m_merge", (
                                 {"type": "multiply"},
                                 {"type": "add"},
                                 {"type": "concatenate", "hidden_fc": None})),
                             },
              "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(0.005)),
              },
    "fit": {"epochs": 150,
            "patience": 2,
            "batch_size": 128,
            }
}

# test_fn(fn, hyper_params, n_train=1000)
trials = CMongoTrials(DB_NAME, ip=host, exp_name, kill_timeout=30 * 60)
best = fmin(fn, hyper_params, trials=trials, algo=tpe.suggest, max_evals=200)
print("best_parameters: " + str(best))
# ----------------
exp_name = "model_deep_position=relu_scaled"
print_exp(exp_name)
# -----
fn = CompileFN(DB_NAME, exp_name,
               data_fn=data.data,
               model_fn=model.model,
               add_eval_metrics=["auprc", "auc"],
               loss_metric="auprc",
               loss_metric_mode="max",
               use_tensorboard=False,
               valid_split=0.2,
               random_state=100)

hyper_params = {
    "data": {"pos_class_weight": 2.0,
             "truncate": False,  # we haven't used truncate either in model_deep_position
             "encode_splines": False,
             "minmax_scale": True,
             },
    "model": {"nonlinearity": "relu",
              "filters": 15 * 2**hp.randint("filters", 5),
              "init_motifs": None,
              "pos_effect": {"type": "relu",
                             "n_bases": 10 + 10 * hp.randint("n_bases", 3),
                             "l2": hp.loguniform("m_pos_effect_l2",
                                                 np.log(1e-18),
                                                 np.log(1e2)),
                             "activation": hp.choice("m_pos_activation", (None, "relu")),
                             "use_bias": False,
                             "merge": {
                                 "type": "concatenate",
                                 "hidden_fc": {"n_hidden": 10 + 30 * hp.randint("n_hidden", 3),
                                               "activation": "relu",
                                               "dropout_rate": hp.uniform("m_dropout_rate", 0, 0.6),
                                               "n_layers": hp.randint("m_n_layers", 4),
                                               },
                             },
                             },
              "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(0.005)),
              },
    "fit": {"epochs": 150,
            "patience": 4,
            "batch_size": 128,
            }
}

# test_fn(fn, hyper_params, n_train=1000)
trials = CMongoTrials(DB_NAME, ip=host, exp_name, kill_timeout=30 * 60)
best = fmin(fn, hyper_params, trials=trials, algo=tpe.suggest, max_evals=1000)
print("best_parameters: " + str(best))

# ------------------------
exp_name = "model_deep_seqonly"
print_exp(exp_name)
# -----
fn = CompileFN(DB_NAME, exp_name,
               data_fn=data.data,
               model_fn=model.model,
               add_eval_metrics=["auprc", "auc"],
               loss_metric="auprc",
               loss_metric_mode="max",
               use_tensorboard=False,
               valid_split=0.2,
               random_state=100)

hyper_params = {
    "data": {"pos_class_weight": 2.0,
             "truncate": False,  # we haven't used truncate either in model_deep_position
             "encode_splines": False,
             "minmax_scale": True,
             },
    "model": {"nonlinearity": "relu",
              "filters": 15 * 2**hp.randint("filters", 5),
              "init_motifs": None,
              "pos_effect": None,
              "lr": hp.loguniform("m_lr", np.log(1e-4), np.log(0.005)),
              },
    "fit": {"epochs": 150,
            "patience": 4,
            "batch_size": 128,
            }
}

test_fn(fn, hyper_params, n_train=1000)
trials = CMongoTrials(DB_NAME, ip=host, exp_name, kill_timeout=30 * 60)
best = fmin(fn, hyper_params, trials=trials, algo=tpe.suggest, max_evals=1000)
print("best_parameters: " + str(best))

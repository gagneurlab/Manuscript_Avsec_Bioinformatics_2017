"""Train a models using different data subsets"""
import os
import argparse
import pandas as pd
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
# my libs
import concise.eval_metrics as cem
from concise.utils.helper import write_json
from concise.hyopt import (CMongoTrials, _train_and_eval_single,
                           get_model, get_data)
# project code
from train_all import DB_NAME, PROC_DIR, POS_FEATURES, DIR_ROOT
from mongodb_setup import host, port
from helper import get_logger
import data
import model


data_fn = data.data_extended_cached
model_fn = model.model

add_eval_metrics = {"auc": cem.auc,
                    "auprc": cem.auprc}
SAVE_ROOT = os.path.join(PROC_DIR, "feature_exclusion_exp")


def get_path(dir_name, file_format, args):
    """Get a normalized path and create a directory if
    it doesn't exist
    """
    fname = "{exp}-excl-{excl}".format(exp=args.exp,
                                       excl=args.feature_set)
    path = os.path.join(SAVE_ROOT, dir_name, args.rbp,
                        fname + file_format)

    # make the directory if it doesn't exist_ok
    os.makedirs(os.path.dirname(path), exist_ok=True)

    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run a program by excluding features")
    parser.add_argument("--rbp", help="Which rbp to use (say UPF1)")
    parser.add_argument("--feature_set",
                        help="Feature set to exclude, separate by ','")
    # with defaults
    parser.add_argument("--exp", default="DeepNN_scalar_position_ext_gam",
                        help="experiment base name to use")

    args = parser.parse_args()

    # output files
    model_path = get_path("models", ".h5", args)
    results_path = get_path("results", ".json", args)
    test_predictions_path = get_path("test_predictions", ".csv", args)
    log_path = get_path("logs", ".log", args)

    logger = get_logger("excl_f", log_path)

    # features need to be in POS_FEATURES
    provided_features = args.feature_set.split(",")
    for f in provided_features:
        assert f in POS_FEATURES
    excl_features = provided_features
    use_features = [f for f in POS_FEATURES if f not in excl_features]

    logger.info("used_features: {0}".format(use_features))

    logger.info("get the best hyper-parameters for a model")
    c_exp_name = args.exp + "_" + args.rbp
    logger.info("c_exp_name: {0}".format(c_exp_name))

    tr = CMongoTrials(DB_NAME, c_exp_name, ip=host, port=port)
    tid = tr.best_trial_tid()
    best_param = tr.get_param(tid)

    # modify which features are we taking
    best_param["model"]["external_pos"]["feat_names"] = use_features

    logger.info("Load the data")
    train, valid, test = get_data(data_fn, best_param)

    logger.info("Initialize the model")
    m = get_model(model_fn, train, best_param)

    logger.info("Setup callbacks")
    callbacks = [EarlyStopping(monitor=best_param["fit"]["early_stop_monitor"],
                               patience=best_param["fit"]["patience"]),
                 ModelCheckpoint(model_path,
                                 monitor=best_param["fit"][
                                     "early_stop_monitor"],
                                 save_best_only=True)]

    logger.info("Train and eval the model on the validation set")
    eval_valid_metrics, history = _train_and_eval_single(train, valid,
                                                         model=m,
                                                         batch_size=best_param[
                                                             "fit"]["batch_size"],
                                                         epochs=best_param[
                                                             "fit"]["epochs"],
                                                         use_weight=False,
                                                         callbacks=callbacks,
                                                         eval_best=True,
                                                         add_eval_metrics=add_eval_metrics)
    # Test-set eval
    logger.info("Load the trained model again")
    trained_model = load_model(model_path)

    logger.info("Predict for the test set")
    y_pred = m.predict(test[0], verbose=0)
    y_true = test[1]

    logger.info("evaluate the test accuracy")
    eval_test_metrics = {k: v(y_true, y_pred)
                         for k, v in add_eval_metrics.items()}

    logger.info("Save the test-predictions")
    dt_pred = pd.DataFrame({"y_true": y_true.reshape((-1,)),
                            "y_pred": y_pred.reshape((-1,))})
    dt_pred.to_csv(test_predictions_path)

    # Final results
    logger.info("Build a results dict")
    ret = {
        # evaluation info
        "eval_valid": eval_valid_metrics,
        "eval_test": eval_test_metrics,
        # features
        "use_features": use_features,
        "excl_features": excl_features,
        # scripts
        "parsed_args": vars(args),
        # experiment info
        "tid": tid,  # trial-id
        "param": best_param,
        "c_exp_name": c_exp_name,
        "exp": args.exp,
        "rbp": args.rbp,
        # save paths
        "path": {
            "model": model_path,
            "results": results_path,
            "test_predictions": test_predictions_path,
        },
        # training history
        "history": history,
    }
    logger.info("Write the results")
    write_json(ret, results_path)

    logger.info("Done!")

"""
file: ml_training_xor_application.py
brief: script to run basic training XOR application using the hipe4ml package
usage: python3 ml_training_xor_application.py --config CONFIG (--train XOR --apply)
"""

import argparse
import os
import pickle
import sys

import uproot

# pylint: disable=import-error
try:
    from hipe4ml import plot_utils
    from hipe4ml.model_handler import ModelHandler
    from hipe4ml.tree_handler import TreeHandler
except ModuleNotFoundError:
    print("Module 'hipe4ml' is not installed. Please install it to run this macro")

import matplotlib.pyplot as plt  # pylint: disable=import-error
import numpy as np  # pylint: disable=import-error
import pandas as pd  # pylint: disable=import-error
import xgboost as xgb  # pylint: disable=import-error
import yaml  # pylint: disable=import-error
from sklearn.model_selection import \
    train_test_split  # pylint: disable=import-error

try:
    from hipe4ml_converter.h4ml_converter import H4MLConverter
except ModuleNotFoundError:
    print("Module 'hipe4ml_converter' is not installed. Please install it to run this macro")

MAX_BKG_FRAC = 0.4  # max of bkg fraction to keep for training


def enforce_list(x):
    """
    Helper method to enforce list type

    Parameters
    -----------------
    - x: a string or a list of string

    Returns
    -----------------
    - x_list if x was not a list (and not None), x itself otherwise
    """

    if not isinstance(x, list):
        # handle possible spaces in config file entry
        x = x.split(",")
        for i, element in enumerate(x):
            x[i] = element.strip()

    return x


# pylint: disable= too-few-public-methods
class MlCommon:
    """
    Class to define common features needed either for ML training&testing or ML application
    """

    def __init__(self, config):
        """
        Init method

        Parameters
        -----------------
        - config: dictionary with config read from a yaml file
        """

        # channel and labels
        self.channel = config["channel"]
        self.labels = enforce_list(config["labels"])
        # input
        self.tree_name = config["tree_name"]
        self.folder_name = config["folder_name"]
        # pt binning
        pt_bins_limits = enforce_list(config["pt_bins_limits"])
        self.pt_bins = [[a, b] for a, b in zip(pt_bins_limits[:-1], pt_bins_limits[1:])]
        # keeping mass and pT
        self.name_pt_var = config["name_pt_var"]
        self.column_to_save_list = config["column_to_save_list"]


# pylint: disable= too-many-instance-attributes, too-few-public-methods
class MlTraining(MlCommon):
    """
    Class for ML training and testing
    """

    def __init__(self, config):
        """
        Init method

        Parameters
        -----------------
        - config: dictionary with config read from a yaml file
        """

        config_common, config_train = config["common"], config["train_ml"]

        # initialize mother common class
        super().__init__(config_common)

        # input
        if config_train["pt_bins_limits"] is not None:
            pt_bins_limits = enforce_list(config_train["pt_bins_limits"])
            self.pt_bins = [[a, b] for a, b in zip(pt_bins_limits[:-1], pt_bins_limits[1:])]
        self.infile_names = config_train["input"]
        self.queries = config_train["query"]

        self.seed_split = config_train["seed_split"]
        # (hyper)parameters
        self.share = config_train["class_balance"]["share"]
        self.bkg_factor = config_train["class_balance"]["bkg_factor"]
        self.downsample_bkg_factor = config_train["downsample_bkg_factor"]
        self.training_vars = enforce_list(config_train["training"]["training_vars"])
        self.test_frac = config_train["test_fraction"]
        self.vars_to_draw = config_train["plots"]["extra_columns"] + self.training_vars

        self.raw_output = config_train["training"]["raw_output"]
        self.roc_auc_approach = config_train["training"]["roc_auc_approach"]
        self.roc_auc_average = config_train["training"]["roc_auc_average"]
        self.score_metric = "roc_auc"
        self.hyper_pars = config_train["training"]["hyper_pars"]
        self.hyper_pars_opt = config_train["training"]["hyper_pars_opt"]

        # output
        self.outdir = config_train["output"]["dir"]
        self.extension = config_train["plots"]["extension"]
        self.log_file = config_train["output"]["log_file"]

    def __check_input_consistency(self):
        """
        Helper method to check self consistency of inputs
        """

        # inputs and queries, defined per class
        for label in self.labels:
            if label not in self.infile_names:
                print(f"\033[91mERROR: no input file defined for the '{label}' class!\033[0m")
                sys.exit()
            if label not in self.queries:
                print(f"\033[91mERROR: no query defined for the '{label}' class!\033[0m")
                sys.exit()
        if "common" not in self.queries:
            print("\033[91mERROR: the 'common' query must be defined (set it to null if not needed)!\033[0m")
            sys.exit()
        # class balance
        if self.share not in ("equal", "all_signal"):
            print(f"\033[91mERROR: class_balance option {self.share} not implemented\033[0m")
            sys.exit()
        if self.share == "all_signal" and len(self.bkg_factor) != len(self.pt_bins):
            print("\033[91mERROR: bkg_factor must be defined for each pT bin!\033[0m")
            sys.exit()
        # training
        if self.training_vars is None:
            print("\033[91mERROR: training columns must be defined!\033[0m")
            sys.exit()
        # hyper-parameters options
        if not isinstance(self.hyper_pars, list):
            print("\033[91mERROR: hyper-parameters must be defined or be a list containing an empty dict!\033[0m")
            sys.exit()
        if not isinstance(self.hyper_pars[0], dict):
            print("\033[91mERROR: hyper-parameters must be a list of dict!\033[0m")
            sys.exit()
        if len(self.hyper_pars) != len(self.pt_bins):
            print("\033[91mERROR: hyper-parameters definition does not match pT binning!\033[0m")
            sys.exit()
        if not isinstance(self.hyper_pars_opt["hyper_par_ranges"], dict):
            print("\033[91mERROR: hyper_pars_opt_config must be defined!\033[0m")
            sys.exit()

    def __get_sliced_dfs(self):
        """
        Helper method to get pT-sliced dataframes for each class

        Returns
        -----------------
        - handlers: dictionary containing TreeHandler objects for each class
        """

        print("Loading and preparing data files: ...", end="\r")

        handlers = {}
        for label in self.labels:
            handlers[label] = TreeHandler(
                file_name=self.infile_names[label],
                tree_name=self.tree_name,
                folder_name=self.folder_name
            ).get_subset(
                f"({self.queries[label]}) and ({self.queries['common']})"
                if self.queries["common"] else self.queries[label]
            )
            handlers[label].slice_data_frame(self.name_pt_var, self.pt_bins, True)

        print("Loading and preparing data files: Done!")
        return handlers

    # pylint: disable=too-many-statements, too-many-branches, too-many-arguments, too-many-locals, too-many-statements
    def __data_prep(self, dfs, pt_bin, out_dir, bkg_factor):
        """
        Helper method for pt-dependent data preparation

        Parameters
        -----------------
        - dfs: dictionary containing pandas dataframes for each class
        - pt_bin: pT bin
        - out_dir: output directory
        - bkg_factor: multiplier for n_sig used to determine n_cand_bkg in the 'all_signal' option

        Returns
        -----------------
        - train_test_data: list containing train/test sets and the associated model predictions
        """
        n_cands = {label: len(dfs[label]) for label in self.labels}
        log_available_cands = (
            f"\nNumber of available candidates "
            f"in {pt_bin[0]} < pT < {pt_bin[1]} GeV/c: \n   "
            + "\n   ".join([f"{lab}: {n_cand}" for lab, n_cand in n_cands.items()])
        )

        print(log_available_cands)

        if self.share == "equal":
            n_cand_min = min(n_cands.values())
            bkg_fraction = n_cand_min / n_cands[self.labels[0]]  # bkg fraction to keep for training
            for lab in self.labels:
                n_cands[lab] = n_cand_min
            log_share = (
                "\nKeep the same number of candidates for each class, "
                "chosen as the minimal number of candidates among all classes."
            )

        elif self.share == "all_signal":
            n_cand_bkg = int(min([n_cands[self.labels[0]], sum(n_cands[lab] for lab in self.labels[1:]) * bkg_factor]))
            log_share = (
                f"\nKeep all signal and use {n_cand_bkg} bkg candidates "
                f"for training and testing ({1 - self.test_frac}-{self.test_frac})"
            )
            bkg_fraction = n_cand_bkg / n_cands[self.labels[0]]
            n_cands[self.labels[0]] = n_cand_bkg

        else:
            print(f"\033[91mERROR: class_balance option {self.share} not implemented\033[0m")
            sys.exit()

        print(log_share)

        log_bkg_fraction = (
            "\nFraction of original (i.e. from original dataset) bkg candidates used for ML: "
            f"{100*bkg_fraction*self.downsample_bkg_factor:.2f}%"
        )
        if (1 - self.test_frac) * bkg_fraction * self.downsample_bkg_factor > MAX_BKG_FRAC:
            log_bkg_fraction += (
                f"\n\033[93m\nWARNING: using more than {100*MAX_BKG_FRAC:.0f}% "
                "of original (i.e. from original dataset) bkg available for training!\033[0m"
            )
        print(log_bkg_fraction)

        log_training_cands = (
            "\nNumber of candidates used for training and testing: \n   "
            + "\n   ".join([f"{lab}: {n_cand}" for lab, n_cand in n_cands.items()])
        )

        print(log_training_cands)

        # write logs in log file
        with open(os.path.join(out_dir, self.log_file), "w", encoding="utf-8") as file:
            file.write(log_available_cands)
            file.write(log_share)
            file.write(log_bkg_fraction)
            file.write(log_training_cands)

        df_tot = pd.concat([dfs[label][:n_cands[label]] for label in self.labels], sort=True)

        labels_array = np.concatenate([[i] * n_cands[self.labels[i]] for i in range(len(self.labels))])
        if 0 < self.test_frac < 1:
            train_set, test_set, y_train, y_test = train_test_split(
                df_tot, labels_array, test_size=self.test_frac, random_state=self.seed_split
            )
        else:
            print("ERROR: test_fraction must belong to ]0,1[")
            sys.exit(0)

        train_test_data = [train_set, y_train, test_set, y_test]
        del df_tot  # release memory

        # safety
        if len(np.unique(train_test_data[3])) != len(self.labels):
            print(
                "\033[91mERROR: The number of labels defined does not match"
                "the number of classes! \nCheck the CONFIG file\033[0m"
            )
            sys.exit()

        # _____________________________________________
        plot_utils.plot_distr(
            [df for df in dfs.values()], self.vars_to_draw, 100, self.labels, figsize=(12, 7), alpha=0.3, log=True, grid=False, density=True
        )
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        for ext in self.extension:
            plt.savefig(f"{out_dir}/DistributionsAll_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        plt.close("all")
        # _____________________________________________
        corr_matrix_fig = plot_utils.plot_corr([df for df in dfs.values()], self.vars_to_draw, self.labels)
        for fig, lab in zip(corr_matrix_fig, self.labels):
            plt.figure(fig.number)
            plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
            for ext in self.extension:
                fig.savefig(f"{out_dir}/CorrMatrix{lab}_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")

        return train_test_data

    # pylint: disable=too-many-statements, too-many-branches
    def __train_test(self, train_test_data, hyper_pars, pt_bin, out_dir):
        """
        Helper method for model training and testing

        Parameters
        -----------------
        - train_test_data: list containing train/test sets and the associated model predictions
        - hyper_pars: default hyper-parameters (can be modified if Optuna enabled)
        - pt_bin: pT bin
        - out_dir: output directory
        """

        model_clf = xgb.XGBClassifier(use_label_encoder=False)
        model_hdl = ModelHandler(model_clf, self.training_vars, hyper_pars)

        # hyperparams optimization
        if self.hyper_pars_opt["activate"]:
            print("Performing optuna hyper-parameters optimisation: ...", end="\r")

            with open(os.path.join(out_dir, self.log_file), "a", encoding="utf-8") as file:
                file.write("\nOptuna hyper-parameters optimisation:")
                sys.stdout = file
                cross_val_scoring = self.score_metric
                if len(self.labels) > 2:
                    cross_val_scoring += f"_{self.roc_auc_approach}"
                    if self.roc_auc_average == "weighted":
                        cross_val_scoring += f"_{self.roc_auc_average}"
                model_hdl.optimize_params_optuna(
                    train_test_data,
                    self.hyper_pars_opt["hyper_par_ranges"],
                    cross_val_scoring=cross_val_scoring,
                    timeout=self.hyper_pars_opt["timeout"],
                    n_jobs=self.hyper_pars_opt["njobs"],
                    n_trials=self.hyper_pars_opt["ntrials"],
                    direction="maximize",
                )
            sys.stdout = sys.__stdout__
            print("Performing optuna hyper-parameters optimisation: Done!")
            print(f"Optuna hyper-parameters:\n{model_hdl.get_model_params()}")
        else:
            model_hdl.set_model_params(hyper_pars)

        # store final hyperparameters in info file
        with open(os.path.join(out_dir, self.log_file), "a", encoding="utf-8") as file:
            file.write(f"\nModel hyperparameters:\n {model_hdl.get_model_params()}")

        # train and test the model with the updated hyper-parameters
        y_pred_test = model_hdl.train_test_model(
            train_test_data,
            True,
            output_margin=self.raw_output,
            average=self.roc_auc_average,
            multi_class_opt=self.roc_auc_approach
        )

        y_pred_train = model_hdl.predict(train_test_data[0], self.raw_output)

        # Save applied model to test set
        test_set_df = train_test_data[2]
        test_set_df = test_set_df.loc[:, self.column_to_save_list]
        if y_pred_test.ndim == 2:
            # For multi-class, store each class probability as a separate column
            for label in self.labels:
                test_set_df[f"ML_output_{label}"] = y_pred_test[:, self.labels.index(label)]
        else:
            # For binary classification
            test_set_df["ML_output"] = y_pred_test

        test_set_df["Labels"] = train_test_data[3]

        if len(self.labels) <= 2:
            suffixes = ["bkg", "signal"]
        else:
            suffixes = self.labels
        for i_label, suffix in enumerate(suffixes):
            test_set_df[test_set_df["Labels"] == i_label].to_parquet(
                f"{out_dir}/{self.channel}_ModelApplied_pT_{pt_bin[0]}_{pt_bin[1]}_{suffix}.parquet.gzip"
            )

        # save model
        if os.path.isfile(f"{out_dir}/ModelHandler_{self.channel}.pickle"):
            os.remove(f"{out_dir}/ModelHandler_{self.channel}.pickle")
        if os.path.isfile(f"{out_dir}/ModelHandler_onnx_{self.channel}.onnx"):
            os.remove(f"{out_dir}/ModelHandler_onnx_{self.channel}.onnx")

        model_hdl.dump_model_handler(f"{out_dir}/ModelHandler_{self.channel}" f"_pT_{pt_bin[0]}_{pt_bin[1]}.pickle")
        model_conv = H4MLConverter(model_hdl)
        model_conv.convert_model_onnx(1)
        model_conv.dump_model_onnx(f"{out_dir}/ModelHandler_onnx_{self.channel}" f"_pT_{pt_bin[0]}_{pt_bin[1]}.onnx")

        # plots
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 7)
        fig_ml_output = plot_utils.plot_output_train_test(
            model_hdl, train_test_data, 80, self.raw_output, self.labels, True, density=True
        )
        for ext in self.extension:
            if isinstance(fig_ml_output, list):
                for fig, lab in zip(fig_ml_output, self.labels):
                    fig.savefig(f"{out_dir}/MLOutputDistr_{lab}_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
            else:
                fig_ml_output.savefig(f"{out_dir}/MLOutputDistr_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 9)
        fig_roc_curve = plot_utils.plot_roc(
            train_test_data[3], y_pred_test, None, self.labels, self.roc_auc_average, self.roc_auc_approach
        )
        for ext in self.extension:
            fig_roc_curve.savefig(f"{out_dir}/ROCCurveAll_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        with open(f"{out_dir}/ROCCurveAll_pT_{pt_bin[0]}_{pt_bin[1]}.pkl", "wb") as file:
            pickle.dump(fig_roc_curve, file)
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 9)
        fig_roc_curve_tt = plot_utils.plot_roc_train_test(
            train_test_data[3],
            y_pred_test,
            train_test_data[1],
            y_pred_train,
            None,
            self.labels,
            self.roc_auc_average,
            self.roc_auc_approach,
        )

        fig_roc_curve_tt.savefig(f"{out_dir}/ROCCurveTrainTest_pT_{pt_bin[0]}_{pt_bin[1]}.pdf")
        # _____________________________________________
        precision_recall_fig = plot_utils.plot_precision_recall(train_test_data[3], y_pred_test, self.labels)
        precision_recall_fig.savefig(f"{out_dir}/PrecisionRecallAll_pT_{pt_bin[0]}_{pt_bin[1]}.pdf")
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (12, 7)
        fig_feat_importance = plot_utils.plot_feature_imp(
            train_test_data[2][train_test_data[0].columns], train_test_data[3], model_hdl, self.labels
        )
        if len(self.labels) <= 2:
            fig_names = [f"FeatureImportance_{self.channel}"]
        else:
            fig_names = [f"FeatureImportance_{lab}_{self.channel}" for lab in self.labels]
        fig_names.append(f"FeatureImportanceAll_{self.channel}")
        for fig, fig_name in zip(fig_feat_importance, fig_names):
            for ext in self.extension:
                fig.savefig(f"{out_dir}/{fig_name}.{ext}")

    def process(self):
        """
        Process function of the class, performing data preparation,
        training, testing, saving the model and important plots
        """

        self.__check_input_consistency()
        dfs = self.__get_sliced_dfs()

        for i_pt, pt_bin in enumerate(self.pt_bins):
            print(f"\n\033[94mStarting ML analysis --- {pt_bin[0]} < pT < {pt_bin[1]} GeV/c\033[0m")

            out_dir_pt = os.path.join(os.path.expanduser(self.outdir), f"pt{pt_bin[0]}_{pt_bin[1]}")
            if os.path.isdir(out_dir_pt):
                print(
                    (
                        f"\033[93mWARNING: Output directory '{out_dir_pt}' already exists,"
                        " overwrites possibly ongoing!\033[0m"
                    )
                )
            else:
                os.makedirs(out_dir_pt)

            if self.share == "all_signal":
                bkg_factor = self.bkg_factor[i_pt]
            else:
                bkg_factor = None

            sliced_dfs = {label: df.get_slice(i_pt) for label, df in dfs.items()}
            train_test_data = self.__data_prep(
                sliced_dfs, pt_bin, out_dir_pt, bkg_factor
            )
            self.__train_test(train_test_data, self.hyper_pars[i_pt], pt_bin, out_dir_pt)


# pylint: disable= too-many-instance-attributes, too-few-public-methods
class MlApplication(MlCommon):
    """
    Class for ML models application
    """

    def __init__(self, config):
        """
        Init method

        Parameters
        -----------------
        - config: dictionary with config read from a yaml file
        """

        config_common, config_apply = config["common"], config["apply_ml"]

        # initialize mother common class
        super().__init__(config_common)

        self.infile_names = config_apply["input"]["file_names"]
        self.model_names = enforce_list(config_apply["input"]["model_names"])
        self.merge_mc_with_check_decay = config_apply["input"]["merge_mc_with_check_decay"]
        self.tree_name_check_decay = config_apply["input"]["tree_name_check_decay"]
        self.outdir = config_apply["output"]["dir"]
        self.out_tree_name = config_apply["output"]["tree_name"]
        self.data_tags = config_apply["output"]["data_tags"]

    def __check_input_consistency(self):
        """
        Helper method to check self consistency of inputs
        """

        if len(self.pt_bins) != len(self.model_names):
            print("\033[91mERROR: pT binning does not match the number of BDT models!\033[0m")
            sys.exit()

    def __load_models(self):
        """
        Helper method to load models

        Returns
        -----------------
        - model_hdls: list of ModelHanlder instances
        """
        model_hdls = []
        for i_bin, _ in enumerate(self.pt_bins):
            path_model = self.model_names[i_bin]
            if not isinstance(path_model, str):
                print("\033[91mERROR: path to model not correctly defined!\033[0m")
                sys.exit()
            path_model = os.path.expanduser(path_model)
            print(f"Loaded saved model: {path_model}")
            model_hdl = ModelHandler()
            model_hdl.load_model_handler(path_model)
            model_hdls.append(model_hdl)
        return model_hdls

    def process(self):
        """
        Process function
        """

        self.__check_input_consistency()
        model_hdls = self.__load_models()

        for infile_name, data_tag in zip(self.infile_names, self.data_tags):
            print(f"Loading and preparing data file {infile_name}: ...", end="\r")
            hdl_data = TreeHandler(file_name=infile_name, tree_name=self.tree_name, folder_name=self.folder_name)
            cols_to_merge = []
            if self.merge_mc_with_check_decay:
                try:
                    hdl_data_check_decay = TreeHandler(
                        file_name=infile_name, tree_name=self.tree_name_check_decay, folder_name=self.folder_name
                    )
                    cols_to_merge = ["fPdgCodeBeautyMother", "fPdgCodeCharmMother"]
                    hdl_data.set_data_frame(pd.concat(
                        [hdl_data.get_data_frame(), hdl_data_check_decay.get_data_frame()[cols_to_merge]],
                        axis=1
                    ))
                except:
                    print(
                        "No deacy check tree found, only the main tree will be used for the application"
                    )
                    cols_to_merge = []
            hdl_data.slice_data_frame(self.name_pt_var, self.pt_bins, True)
            print(f"Loading and preparing data files {infile_name}: Done!")

            out_dir = os.path.expanduser(self.outdir)
            if os.path.isdir(out_dir):
                print(
                    (
                        f"\033[93mWARNING: Output directory '{out_dir}' already exists,"
                        " overwrites possibly ongoing!\033[0m"
                    )
                )
            else:
                os.makedirs(out_dir)
            print("Applying ML model to dataframes: ...", end="\r")
            for ibin, pt_bin in enumerate(self.pt_bins):
                df_data_pt_sel = hdl_data.get_slice(ibin)
                ypred = model_hdls[ibin].predict(df_data_pt_sel, False)

                df_data_pt_sel = df_data_pt_sel.loc[:, self.column_to_save_list + cols_to_merge]
                if ypred.ndim == 2:
                    # Multi-class
                    for i_label, label in enumerate(self.labels):
                        df_data_pt_sel[f"ML_output_{label}"] = ypred[:, i_label]
                else:
                    # Binary
                    df_data_pt_sel["ML_output"] = ypred

                outfile_name = f"{out_dir}/{data_tag}_{self.channel}_pT_{pt_bin[0]}_{pt_bin[1]}_ModelApplied"
                outfile_name_root = outfile_name + ".root"
                with uproot.recreate(outfile_name_root) as ofile:
                    ofile[self.out_tree_name] = df_data_pt_sel
                outfile_name_parquet = outfile_name + ".parquet.gzip"
                df_data_pt_sel.to_parquet(outfile_name_parquet)

                del df_data_pt_sel


def main(cfg, train):
    """
    Main function

    Parameters
    -----------------
    - cfg: dictionary with config read from a yaml file

    - train: boolean implying ML training&testing (if true) or application (if false)
    """

    if train:
        MlTraining(cfg).process()
    else:  # if we do not train, we apply
        MlApplication(cfg).process()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-c", metavar="text", default="config_ml.yml",
                        help="yaml config file for BDT training and application", required=True)
    parser.add_argument("--train", help="perform only training and testing", action="store_true")
    parser.add_argument("--apply", help="perform only application", action="store_true")
    args = parser.parse_args()

    print("Loading configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading configuration: Done!")

    if args.train and args.apply:
        print("\033[91mERROR: Both --train and --apply options are activated, choose only one.\033[0m")
        sys.exit()
    elif not args.train and not args.apply:
        print("\033[91mERROR: None of --train or --apply options are activated, choose (only) one.\033[0m")
        sys.exit()

    main(configuration, args.train)

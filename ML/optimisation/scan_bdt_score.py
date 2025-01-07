"""
Script to scan B-meson BDT score and compute relevant quantities to optimise selections
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot
import yaml
from hist import Hist

sys.path.append('../../utils')
os.environ["CUDA_VISIBLE_DEVICES"] = ""
from df_utils import read_parquet_in_batches
from analysis_utils import get_n_events_from_zorro
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
import zfit

# TODO currently, the script assumes that the input distributions have the same pT binning as the one defined in
# the configuration file !

# FIXME some utils can't be taken from utils file if PyROOT has compatibility issues
# maybe add a try except test in the utils files like shown below:
# try:
#     import ROOT
# except ImportError:
#     print("\033[93mWARNING: ROOT is not installed or has not been built for this version of Python!\033[0m")
#     print("\033[93mIgnore this message if no ROOT feature is used in the running script.\033[0m")


def enforce_trailing_slash(path):
    """
    Helper method to enforce '/' at the and of directory name

    Parameters
    ----------
    - path: str
        Some path

    Returns
    ----------
    - path: str
        Path with a trailing slash at the end if it was not there yet
    """

    if path is not None and path[-1] != "/":
        path += "/"

    return path


# FIXME this method is in the extract_raw_yield macro, we should move it to Utils/
def create_hist(x_axis, contents, errors, label="BDT score"):
    """
    Helper method to create histogram

    Parameters
    ----------
    - x_axis (list): x axis list
    - contents (list): histogram contents
    - errors (list): histogram errors
    - label (str): label for x axis

    Returns
    ----------
    - histogram (hist.Hist)

    """
    x_bin_cent = [0.5*(bdt_min + bdt_max) for bdt_min, bdt_max in zip(x_axis[:-1], x_axis[1:])]
    histo = Hist.new.Var(x_axis, name="x", label=label).Weight()
    histo.fill(x_bin_cent, weight=contents)
    histo.view(flow=False).variance = np.array(errors)**2

    return histo


def get_labeled_dfs(infile_names):
    """
    Helper method to get labeled dataframes from input files

    Parameters
    ----------
    - infile_names: list(str)
        Names of input files

    Returns
    ----------
    - df_mc_prd_bkg: pandas.DataFrame
        Partly reconstructed decays background-labeled dataframe

    - df_mc_sig: pandas.DataFrame
        Signal-labeled dataframe
    """

    # load dataframes from input files
    df_tot = pd.concat([read_parquet_in_batches(parquet) for parquet in infile_names])
    df_mc_prd_bkg = df_tot.query("fFlagMcMatchRec == 8")
    df_mc_sig = df_tot.query("fFlagMcMatchRec == -1 or fFlagMcMatchRec == 1")

    return df_mc_prd_bkg, df_mc_sig


# TODO add possibility to import FONLL graph with different binning
def get_fonll_cross_section(file_name, graph_name):
    """
    Helper method to get FONLL cross section

    Parameters
    -----------------
    - file_name: str
        Name of the input file with FONLL information

    - graph_name: str
        Name of the graph in the input file

    Returns
    -----------------
    - values: list(float)
        FONLL pT-differential cross section values

    - low_uncs: list(float)
        Lower uncertainties on FONLL pT-differential cross section values

    - high_uncs: list(float)
        Higher uncertainties on FONLL pT-differential cross section values
    """
    with uproot.open(file_name)[graph_name] as graph:
        values = graph.values(axis="y")
        low_uncs = graph.errors(which="low", axis="y")
        high_uncs = graph.errors(which="high", axis="y")

    return values, low_uncs, high_uncs


def get_integrated_luminosity(
              infile_names,
              zorro_folder,
              triggers_of_interest_names,
              h_collisions_path,
              tvx_cross_section
        ):
    """
    Helper method to get the integrated luminosity

    Parameters
    -----------------
    - name_analysis_results_files: str or list(str)
        Name(s) of the .root input file(s)

    - name_zorro_folder: str
        Name of Zorro folder

    - triggers_of_interest_names: str or list(str)
        Name(s) of the triggers of interest

    - h_collisions_path: str
        Path to the hCollisions histogram

    - tvx_cross_section: float
        TVX cross section

    Returns
    -----------------
    - int_lumi: float
        Integrated luminosity
    """
    n_events = get_n_events_from_zorro(
        infile_names,
        zorro_folder,
        triggers_of_interest_names,
        h_collisions_path
        )
    int_lumi = float(n_events) / tvx_cross_section
    return int_lumi


def get_expected_signal_over_deltapt_over_acceff(fonll_cross_section_tuple, frag_frac, br, int_lumi):
    """
    Helper method to retrieve the expected signal divided by pT bin width and (Acc x eff)

    Parameters
    -----------------
    - fonll_cross_section_tuple: (list(float), list(float), list(float))
        Tuple containing FONLL predictions on B hadron cross section (values and uncertainties)

    - frag_frac: float
        Fragmentation fraction (B hadron -> B0/B+)

    - br: float
        Branching ratio of the decay channel

    - int_lumi: float
        Integrated luminosity

    Returns
    -----------------
    - exp_sig_over_dpt_acceff_tuple: list(list(float, float, float))
        Expected signal divided by pT bin width and (Acc x eff)
        exp_sig_over_dpt_acceff_tuple[ipt] contains the tuple (value, low_unc, high_unc) for a pT bin of index ipt
    """

    coef = 2 * frag_frac * br * int_lumi
    exp_sig_over_dpt_acceff_tuple = []
    for value, low_unc, high_unc in zip(fonll_cross_section_tuple[0],
                                        fonll_cross_section_tuple[1],
                                        fonll_cross_section_tuple[2]):
        exp_sig_over_dpt_acceff_tuple.append([value * coef, low_unc * coef, high_unc * coef])

    return exp_sig_over_dpt_acceff_tuple


def get_expected_signal(exp_sig_over_dpt_acceff_tuple, dpt, acc_eff, acc_eff_unc):
    """
    Helper method to retrieve the expected signal in given pT bin and a given BDT score cut

    Parameters
    -----------------
    - exp_sig_over_dpt_acceff_tuple: list(float, float, float)
        Expected signal divided by pT bin width and (Acc x eff) for a given pT bin

    - dpt: float
        pT bin width

    - acc_eff: float
        Acceptance times efficiency
    
    - acc_eff_unc: float
        Uncertainty on acceptance times efficiency

    - asymm_unc: bool
        Switch to keep asymmetric uncertainties coming from FONLL predictions

    Returns
    -----------------
    - exp_sig: float
        Expected signal divided by pT bin width and (Acc x eff)
    """

    coef = dpt * acc_eff
    exp_sig_value = exp_sig_over_dpt_acceff_tuple[0] * coef
    exp_sig_low_unc = exp_sig_over_dpt_acceff_tuple[1] * coef
    exp_sig_high_unc = exp_sig_over_dpt_acceff_tuple[2] * coef

    exp_sig_unc = (exp_sig_low_unc + exp_sig_high_unc) / 2 / exp_sig_value # relative uncertainty
    exp_sig_unc = np.sqrt(exp_sig_unc**2 + (acc_eff_unc / acc_eff)**2) * exp_sig_value

    return exp_sig_value, exp_sig_unc


# TODO add safety by checking edges for pT binning compatibility
def get_acc_eff_presel(file_name, hist_name):
    """
    Helper method to retrieve the preselections acceptance times efficiency

    Parameters
    -----------------
    - file_name: str
        Name of the .root file containing preselections acceptance times efficiency

    - hist_name: str
        Name of the TH1 histogram containing preselections acceptance times efficiency

    Returns
    -----------------
    - values, uncs: list(float), list(float)
        Values and uncertainties of preselections acceptance times efficiency
    """
    with uproot.open(file_name)[hist_name] as hist:
        values = hist.values()
        uncs = hist.errors()
        bins = hist.axis().edges()

    return values, uncs, bins


def get_signal_region(df_mc_sig, fit_limits, nbins, pt_bin, sel, pdg_code, verbosity=None):
    """
    Helper method to get B(3sigma) from sidebands

    Parameters
    -----------------
    - df_mc_sig: pandas.DataFrame
        Dataframe containing MC signal

    - fit_limits: list(float)
        Fitting range

    - nbins: int
        Number of bins chosen by user to bin data in case of unbinned data

    - pt_bin: list(float)
        pT bin

    - sel: float
        Cut on the BDT output score

    - pdg_code: int
        PDG code of the beauty meson

    - verbosity: int
        Verbosity level (from 0 to 10) Default value to 0

    Returns
    -----------------
    - signal_region: list(float)
        The region of the signal peak [mu - 3sigma, mu + 3sigma]
    """
    data_hdl_mc = DataHandler(df_mc_sig,
                              var_name="fM",
                              limits=fit_limits,
                              nbins=nbins)
    fitter = F2MassFitter(data_hdl_mc,
                          ["gaussian"],  # we set gaussian as we want to retrieve sigma
                          ["nobkg"],
                          name=f"b0_mc_pT_{pt_bin[0]:.0f}_{pt_bin[1]:.0f}_ML_cut_{1000*sel:.0f}",
                          verbosity=verbosity)
    fitter.set_signal_initpar(0, "sigma", 0.03, limits=[0.01, 0.08])
    fitter.set_particle_mass(0, pdg_id=pdg_code)
    result = fitter.mass_zfit()
    if result.converged:
        mu, _ = fitter.get_mass(0)
        sigma, _ = fitter.get_sigma(0)

    signal_region = [mu - 3 * sigma, mu + 3 * sigma]
    return signal_region


# pylint: disable=too-many-locals,too-many-arguments
def get_expected_background_from_sidebands(
        df_data,
        full_inv_mass_region,
        sidebands_regions,
        signal_region,
        nbins,
        bkg_funcs,
        pt_bin,
        sel,
        outdir_name,
        df_mc_prd_bkg=None,
        verbosity=None):
    """
    Helper method to get B(3sigma) from sidebands

    Parameters
    -----------------
    - df_data: pandas.DataFrame
        Dataframe containing MC signal

    - full_inv_mass_region: list(float)
        The full fitting range considered (sidebands + excluded region under the expected signal peak)

    - sidebands_region: list(list((float))
        The fitting range for sidebands (leftband + rightband)

    - signal_region: list(float)
        Signal region determined by a fit of MC signal

    - nbins: int
        Number of bins chosen by user to bin data in case of unbinned data

    - bkg_funcs: list(str)
        Names of bkg pdfs

    - pt_bin: list(float)
        pT bin

    - sel: float
        Cut on the BDT output score

    - outdir_name: str
        Name of the output directory

    - df_mc_prd_bkg: pandas.DataFrame
        Dataframe containing MC partly reco decays

    - verbosity: int
        Verbosity level (from 0 to 10) Default value to 0

    Returns
    -----------------
    - exp_bkg, exp_bkg_unc: (float, float)
        Expected background and its uncertainty in the signal region
    """
    data_hdl = DataHandler(df_data,
                           var_name="fM",
                           limits=full_inv_mass_region,
                           nbins=nbins)

    if df_mc_prd_bkg is not None:
        data_hdl_prd_bkg = DataHandler(df_mc_prd_bkg,
                                       var_name="fM",
                                       limits=full_inv_mass_region,
                                       nbins=nbins)
        bkg_funcs.insert(0, "kde_grid")

    fitter = F2MassFitter(data_hdl,
                          ["nosignal"],
                          bkg_funcs,
                          name=f"bkg_pT_{pt_bin[0]:.0f}_{pt_bin[1]:.0f}_ML_cut_{1000*sel:.0f}",
                          limits=sidebands_regions,
                          verbosity=verbosity)

    if df_mc_prd_bkg is not None:
        fitter.set_background_kde(0, data_hdl_prd_bkg)

    fitter.set_background_initpar(0, "frac", 0.2, limits=[0., 1.])
    fitter.set_background_initpar(1, "lam", -1.2, limits=[-10., 10.])
    fitter.set_background_initpar(1, "c1", -0.05, limits=[-0.2, 0.])
    fitter.set_background_initpar(1, "c2", 0.008, limits=[0.000, 0.03])
    result = fitter.mass_zfit()

    exp_bkg, exp_bkg_unc = 0., 0.
    if result.converged:
        exp_bkg, exp_bkg_unc = fitter.get_background(min=signal_region[0], max=signal_region[1])

        fig, _ = fitter.plot_mass_fit(
            style="ATLAS",
            figsize=(8, 8),
            axis_title=r"$M(\mathrm{D^-\pi^+})$ (GeV/$c^2$)",
            show_extra_info=False)
        fig.savefig(f"{outdir_name}sidebands_fit_pT_{pt_bin[0]:.0f}_{pt_bin[1]:.0f}_ML_cut_{1000*sel:.0f}.pdf")

    return exp_bkg, exp_bkg_unc


def get_expected_significance(exp_sig, exp_sig_unc, exp_bkg, exp_bkg_unc):
    """
    Helper method to get B(3sigma) from sidebands

    Parameters
    -----------------
    - exp_sig: float
        Expected signal

    - exp_sig_unc: float
        Uncertainty on expected signal

    - exp_bkg: float
        Expected bkg

    - exp_bkg_unc: float
        Uncertainty on expected bkg

    Returns
    -----------------
    - exp_signif, exp_signif_unc: (float, float)
        Expected significances and their uncertainties
    """
    exp_signif = exp_sig / (exp_sig + exp_bkg)**0.5

    partial_deriv_wrt_sig = (0.5*exp_sig + exp_bkg) / (exp_sig + exp_bkg)**1.5
    partial_deriv_wrt_bkg = (-0.5*exp_sig) / (exp_sig + exp_bkg)**1.5
    exp_signif_unc = (partial_deriv_wrt_sig**2 * exp_sig_unc**2 + partial_deriv_wrt_bkg**2 * exp_bkg_unc**2)**0.5

    return exp_signif, exp_signif_unc

# pylint: disable=too-many-locals,too-many-statements,too-many-branches
def scan(config):
    """
    Main method for BDT score scan

    Parameter
    -----------------
    - config: dict
        Configuration from the YAML file
    """

    # pT bining
    channel = config["channel"]["name"]
    pdg_code = config["channel"]["pdg_code"]
    pt_mins = config['pt']['mins']
    pt_maxs = config['pt']['maxs']
    # configure output
    outdir = enforce_trailing_slash(config["output"]["outdir"])
    sidebands_fit_dir = outdir + channel + "_" + enforce_trailing_slash(config["output"]["sidebands_fit_dir"])
    scan_results_file_name = outdir + channel + "_" + config["output"]["scan_results_file_name"]
    for directory in [outdir, sidebands_fit_dir]:
        if os.path.isdir(directory):
            print(
                (
                    f"\033[93mWARNING: Output directory '{directory}' already exists,"
                    " overwrites possibly ongoing!\033[0m"
                )
            )
        else:
            os.makedirs(directory)

    # compute expected signal over pt bin width and (acceptance times efficiency)
    # TODO add safety to perform pT binning matching
    fonll_cross_section_tuple = get_fonll_cross_section(config["inputs"]["fonll"]["file_name"],
                                                        config["inputs"]["fonll"]["graph_name"])
    # TODO make a safety using the info on tgraph error x axis
    # if len(fonll_cross_section_tuple[0]) != len(pt_mins):
    #     print("\033[91mERROR: FONLL pT binning does not match the chosen pT binning!\033[0m")

    int_lumi = get_integrated_luminosity(config["inputs"]["zorro"]["file_names"],
                                         config["inputs"]["zorro"]["folder_name"],
                                         config["inputs"]["zorro"]["triggers_of_interest"],
                                         config["inputs"]["zorro"]["h_collisions_path"],
                                         config["inputs"]["zorro"]["tvx_cross_section"])

    frag_frac = config["inputs"]["fonll"]["frag_frac"]
    if channel == "B0ToDPi":
        br = config["channel"]["br"]["b0_todminuspi"] * config["channel"]["br"]["dplus_tokpipi"]
    # Get expected signal divided by Delta pT and (Acc x eff)
    # then we only have to multiply by Delta pT and (Acc x eff) to get the expected signal for each BDT score cut
    exp_sig_over_deltapt_over_acceff_tuple = get_expected_signal_over_deltapt_over_acceff(
        fonll_cross_section_tuple,
        frag_frac,
        br,
        int_lumi)

    df_mc_prd_bkg, df_mc_sig = get_labeled_dfs(config["inputs"]["mc"]["file_names"])

    # Retrieve the histogram conatining acceptance times efficiency of preselections
    acc_eff_presel, acc_eff_unc_presel, acc_eff_bins = get_acc_eff_presel(
        config["inputs"]["acc_eff_presel"]["file_name"],
        config["inputs"]["acc_eff_presel"]["hist_name"])

    # Get real data candidates for expected background computation
    df_data = pd.concat(
        [read_parquet_in_batches(parquet) for parquet in config["inputs"]["real_data"]["file_names"]])
    leftband_range = [config["fit_bkg"]["mass_limits_for_fit"][0], config["fit_bkg"]["mass_range_to_exclude"][0]]
    rightband_range = [config["fit_bkg"]["mass_range_to_exclude"][1], config["fit_bkg"]["mass_limits_for_fit"][1]]
    bkg_funcs = config["fit_bkg"]["bkg_funcs"]
    nbins = config["fit_bkg"]["nbins"]
    fit_verbosity = config["fit_bkg"]["verbosity"]

    # loop over pT bins
    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        # configure scan
        selection_steps = config['ML_selections']['steps'][i_pt]
        selection_mins = config['ML_selections']['mins'][i_pt]
        selection_maxs = config['ML_selections']['maxs'][i_pt] + 0.1 * selection_steps
        bdt_selections = np.arange(selection_mins, selection_maxs, selection_steps)
        print(bdt_selections)
        # define output histograms xaxis
        xaxis_hist = (bdt_selections.copy() - selection_steps/2).tolist()
        xaxis_hist.append(bdt_selections[-1] + selection_steps/2)
        bdt_selections = bdt_selections.tolist()
        # get the dataframes entries for this pT interval
        df_mc_sig_pt = df_mc_sig.query(f"{pt_min} < fPt < {pt_max}")
        df_data_pt = df_data.query(f"{pt_min} < fPt < {pt_max}")
        if config["fit_bkg"]["use_bkg_templ"]:
            df_mc_prd_bkg_pt = df_mc_prd_bkg.query(f"{pt_min} < fPt < {pt_max}")

        acc_eff_presel_ipt = acc_eff_presel[np.digitize((pt_min+pt_max)/2, acc_eff_bins) - 1]
        acc_eff_unc_presel_ipt = acc_eff_unc_presel[np.digitize((pt_min+pt_max)/2, acc_eff_bins) - 1]

        # perform the scan
        exp_sig_ipt, exp_sig_unc_ipt = [[] for _ in range(2)]
        exp_bkg_ipt, exp_bkg_unc_ipt = [[] for _ in range(2)]
        acc_eff_ipt, acc_eff_unc_ipt = [[] for _ in range(2)]
        exp_signif_ipt, exp_signif_unc_ipt = [[] for _ in range(2)]
        print(f"Starting ML score scan for {pt_min:.0f} < pT < {pt_max:.0f}: ...", end="\r")
        for i_sel, bdt_sel in enumerate(bdt_selections):
            df_mc_sig_pt_sel = df_mc_sig_pt.query(f"ML_output > {bdt_sel}")
            eff_bdt = len(df_mc_sig_pt_sel) / len(df_mc_sig_pt)  # TODO add uncertainty for BDT efficiency
            acc_eff, acc_eff_unc = acc_eff_presel_ipt * eff_bdt, acc_eff_unc_presel_ipt * eff_bdt

            # expected signal
            exp_sig, exp_sig_unc = get_expected_signal(exp_sig_over_deltapt_over_acceff_tuple[i_pt],
                                                       pt_max-pt_min,
                                                       acc_eff, acc_eff_unc)

            # expected bkg
            signal_region = get_signal_region(df_mc_sig_pt_sel,
                                              config["fit_bkg"]["mass_limits_for_fit"],
                                              nbins,
                                              [pt_min, pt_max],
                                              i_sel,
                                              pdg_code,
                                              verbosity=fit_verbosity)

            if config["fit_bkg"]["use_bkg_templ"]:
                df_for_template_bkg = df_mc_prd_bkg_pt.query(f"ML_output > {bdt_sel}")
                correlated_bkgs = config["fit_bkg"]["correlated_bkgs"]
                dfs_prd_bkg_orig_pt, fracs_pt = [], []
                den_norm = len(df_mc_sig_pt_sel) * config["fit_bkg"]["signal_br"]["pdg"] / config["fit_bkg"]["signal_br"]["sim"]
                for bkg in correlated_bkgs:
                    dfs_prd_bkg_pt = df_for_template_bkg.query("fFlagMcMatchRec == 8 "
                        f"and fPdgCodeBeautyMother == {bkg['beauty_id']} and "
                        f"fPdgCodeCharmMother == {bkg['charm_id']}"
                    )
                    dfs_prd_bkg_orig_pt.append(dfs_prd_bkg_pt)
                    # store all fractions to fix it later in the fit
                    fracs_pt.append(len(dfs_prd_bkg_pt) * bkg["br_pdg"] / bkg["br_sim"] / den_norm)
                # sum of all fractions to normalise the weighted average
                sum_fracs_pt = sum(fracs_pt)
                fracs_pt_norm = [frac / sum_fracs_pt for frac in fracs_pt]

                # dfs for weighted average
                dfs_prd_bkg_sampled = []
                for frac, df_bkg in zip(fracs_pt_norm, dfs_prd_bkg_orig_pt):
                    dfs_prd_bkg_sampled.append(df_bkg.sample(frac=frac))

                df_for_template_bkg_sampled = pd.concat(dfs_prd_bkg_sampled)
            else:
                df_for_template_bkg = None

            exp_bkg, exp_bkg_unc = get_expected_background_from_sidebands(
                df_data_pt.query(f"ML_output > {bdt_sel}"),
                config["fit_bkg"]["mass_limits_for_fit"],
                [leftband_range, rightband_range],
                signal_region,
                nbins,
                bkg_funcs[i_pt].copy(),
                [pt_min, pt_max],
                bdt_sel,
                sidebands_fit_dir,
                df_mc_prd_bkg=df_for_template_bkg_sampled,
                verbosity=fit_verbosity)

            del df_for_template_bkg

            # expected significance
            exp_signif, exp_signif_unc = get_expected_significance(exp_sig,
                                                                   exp_sig_unc,
                                                                   exp_bkg,
                                                                   exp_bkg_unc)

            # fill lists
            exp_sig_ipt.append(exp_sig)
            exp_sig_unc_ipt.append(exp_sig_unc)
            exp_bkg_ipt.append(exp_bkg)
            exp_bkg_unc_ipt.append(exp_bkg_unc)
            acc_eff_ipt.append(acc_eff)
            acc_eff_unc_ipt.append(acc_eff_unc)
            exp_signif_ipt.append(exp_signif)
            exp_signif_unc_ipt.append(exp_signif_unc)

        values_quantities = [exp_sig_ipt, exp_bkg_ipt, acc_eff_ipt, exp_signif_ipt]
        uncs_quantities = [exp_sig_unc_ipt, exp_bkg_unc_ipt, acc_eff_unc_ipt, exp_signif_unc_ipt]
        hists = []
        for values_quantity, uncs_quantity in zip(values_quantities, uncs_quantities):
            hists.append(create_hist(xaxis_hist, values_quantity, uncs_quantity, label="BDT score"))

        ofile_dir_name = f"pT_{pt_min:.0f}_{pt_max:.0f}"
        labels = ["ExpSig", "ExpBkg", "AccEff", "ExpSignif"]
        if i_pt == 0:
            with uproot.recreate(scan_results_file_name) as ofile:
                for label, hist in zip(labels, hists):
                    ofile[f"{ofile_dir_name}/h_{label}_vs_sel"] = hist
        else:
            with uproot.update(scan_results_file_name) as ofile:
                for label, hist in zip(labels, hists):
                    ofile[f"{ofile_dir_name}/h_{label}_vs_sel"] = hist

        # fill canvas with all the relevant quantities
        ylabels = ["expected signal",
                   "expected background",
                   "Acc x eff",
                   "expected significance"]
        fig, axs = plt.subplots(2, 2, figsize=(10, 8), tight_layout=True)
        for ihist, (hist, ylabel) in enumerate(zip(hists, ylabels)):
            ax = axs[0][ihist] if ihist < 2 else axs[1][ihist-2]
            ax.set_ylabel(ylabel)
            ax.set_xlim(xaxis_hist[0], xaxis_hist[-1])
            hist.plot(ax=ax)
        fig.savefig(f"{outdir}scan_{ofile_dir_name}.png")
        print(f"Starting ML score scan for {pt_min:.0f} < pT < {pt_max:.0f}: Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-c", metavar="text", default="config_scan.yml",
                        help="yaml config file for BDT score scan", required=True)
    args = parser.parse_args()

    print("Loading configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading configuration: Done!")
    zfit.run.set_cpus_explicit(intra=configuration["zfit_cpus"]["intra"], inter=configuration["zfit_cpus"]["inter"])
    scan(configuration)

"""
This script performs a multitrial procedure to get the raw yields systematic uncertainty.

Usage:
    python get_ry_systematic.py configFile
"""

import re
import argparse
import itertools
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import yaml
import uproot
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Rectangle
import pandas as pd
import zfit
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from flarefly.utils import Logger


def draw_multitrial(df_multitrial, cfg, pt_min, pt_max, idx_assigned_syst, h_rawy, h_sigma):  # pylint: disable=too-many-locals, too-many-statements # noqa: 501
    """
    Produce a plot with the results of the multitrial procedure.

    Parameters:
    - df_multitrial (DataFrame): DataFrame containing the multitrial data.
    - cfg (dict): Configuration dictionary.
    - pt_min (float): Minimum pt value.
    - pt_max (float): Maximum pt value.
    - idx_assigned_syst (int): Index of the assigned systematic uncertainty.

    Returns:
    - None
    """
    multitrial_cfg = cfg["multitrial"]
    n_trials = len(df_multitrial)
    n_bincounts = len(multitrial_cfg['bincounting_nsigma'])
    x_axis_range = n_trials * (n_bincounts + 1) + 1

    fig, axs = plt.subplots(2, 2, figsize=(20, 15))

    i_pt = np.digitize((pt_min + pt_max) / 2, h_rawy.axis().edges()) - 1
    central_rawy = h_rawy.values()[i_pt]
    central_rawy_unc = h_rawy.errors()[i_pt]
    central_sigma = h_sigma.values()[i_pt]
    central_sigma_unc = h_sigma.errors()[i_pt]

    # Plot the results
    axs[0, 0].errorbar(
        x=range(1, n_trials + 1), y=df_multitrial["rawy"],
        yerr=df_multitrial["rawy_unc"], fmt='o', label='Fit',
        zorder=2
    )

    for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
        axs[0, 0].errorbar(
            x=list(range(n_trials * (i_nsigma + 1) + 1, n_trials * (i_nsigma + 2) + 1)),
            y=df_multitrial[f"rawy_bincounting_{nsigma}"],
            yerr=df_multitrial[f"rawy_bincounting_{nsigma}_unc"], fmt='o',
            label=fr'Bin counting {nsigma}$\sigma$',
            zorder=1
        )

    # Draw the central values
    axs[0, 0].axhline(y=central_rawy, color='r', linestyle='--')
    axs[0, 0].add_patch(
        Rectangle(
            (0, central_rawy - central_rawy_unc),
            x_axis_range, 2 * central_rawy_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'Central value $\pm$ uncertainty'
        )
    )

    axs[0, 0].set_xlim(0, x_axis_range)
    axs[0, 0].set_xlabel('Trial', fontsize=14)
    axs[0, 0].set_ylabel('Raw yield', fontsize=14)
    axs[0, 0].legend(fontsize=12)

    # Draw the raw yields distribution
    axs[0, 1].hist(
        df_multitrial["rawy"], bins=30, alpha=0.7, label='Fit',
        histtype='stepfilled', ec="k", linewidth=2, zorder=2
    )

    for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
        axs[0, 1].hist(
            df_multitrial[f"rawy_bincounting_{nsigma}"],
            bins=30,
            alpha=0.3,
            label=fr'Bin Counting {nsigma}$\sigma$',
            histtype='stepfilled',
            ec="k",
            zorder=1
        )

    # Draw information
    info = 'Fit:\n'
    info += fr'$\mu =$ {np.mean(df_multitrial["rawy"]):.3f}''\n'
    info += fr'$\sigma =$ {np.std(df_multitrial["rawy"]):.3f}''\n'
    for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
        info += fr'{nsigma}$\sigma$ Bin counting:''\n'
        info += fr'$\mu =$ {np.mean(df_multitrial[f"rawy_bincounting_{nsigma}"]):.3f}''\n'
        info += fr'$\sigma =$ {np.std(df_multitrial[f"rawy_bincounting_{nsigma}"]):.3f}''\n'
    anchored_text_fit = AnchoredText(
        info,
        loc='upper left',
        frameon=False
    )

    # Draw the rms + shift from the central value
    rms_shift = get_rms_shift_sum_quadrature(df_multitrial, h_rawy, i_pt)
    axs[0, 1].axvline(x=central_rawy, color='r', linestyle='--')
    axs[0, 1].add_patch(
        Rectangle(
            (central_rawy - rms_shift, 0),
            2 * rms_shift, axs[0, 1].get_ylim()[1],
            color='r', alpha=0.3, zorder=0,
            label=r'$\mathrm{\sqrt{RMS^2 + \Delta^2}}$'
        )
    )
    axs[0, 1].add_artist(anchored_text_fit)

    # Draw the assigned systematic uncertainty
    axs[0, 1].add_patch(
        Rectangle(
            (central_rawy - cfg["assigned_syst"][idx_assigned_syst] * central_rawy, 0),
            2 * cfg["assigned_syst"][idx_assigned_syst] * central_rawy, axs[0, 1].get_ylim()[1],
            color='limegreen', alpha=0.3, zorder=0,
            label='Assigned syst.'
        )
    )

    axs[0, 1].set_xlabel('Raw yields', fontsize=14)
    axs[0, 1].set_ylabel('Counts', fontsize=14)
    axs[0, 1].legend(fontsize=12, loc='upper right')

    x_min = min(
        df_multitrial["rawy"].min(),
        *[df_multitrial[f"rawy_bincounting_{nsigma}"].min()
            for nsigma in multitrial_cfg['bincounting_nsigma']],
        central_rawy - rms_shift
    )
    x_max = max(
        df_multitrial["rawy"].max(),
        *[df_multitrial[f"rawy_bincounting_{nsigma}"].max()
            for nsigma in multitrial_cfg['bincounting_nsigma']],
        central_rawy + rms_shift
    )
    axs[0, 1].set_xlim(x_min * 0.9, x_max * 1.1)

    # Draw the peak width
    axs[1, 0].errorbar(
        x=range(1, len(df_multitrial["sigma"]) + 1),
        y=df_multitrial["sigma"],
        yerr=df_multitrial["sigma_unc"], fmt='o',
        zorder=2
    )

    # Draw the central values
    axs[1, 0].axhline(y=central_sigma, color='r', linestyle='--')
    axs[1, 0].add_patch(
        Rectangle(
            (0, central_sigma - central_sigma_unc),
            x_axis_range, 2 * central_sigma_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'Central value $\pm$ uncertainty'
        )
    )

    axs[1, 0].set_xlim(0, x_axis_range)
    axs[1, 0].set_xlabel('Trial', fontsize=14)
    axs[1, 0].set_ylabel('Width ($GeV/c^2$)', fontsize=14)
    axs[1, 0].legend(fontsize=12)

    # Draw the chi2/ndf
    axs[1, 1].scatter(
        x=range(1, len(df_multitrial["chi2_ndf"]) + 1),
        y=df_multitrial["chi2_ndf"]
    )
    axs[1, 1].set_xlim(0, x_axis_range)
    axs[1, 1].set_xlabel('Trial', fontsize=14)
    axs[1, 1].set_ylabel(r'$\chi^2/$ndf', fontsize=14)

    fig.savefig(
        os.path.join(cfg["output_dir"], f'fig_{pt_min*10:.0f}_{pt_max*10:.0f}.png'),
        bbox_inches='tight'
    )


def get_input_data(cfg, pt_mins, pt_maxs, bdt_cut_mins, bdt_cut_maxs):  # pylint: disable=too-many-locals # noqa: 501
    """
    Retrieve input data based on the given configuration and selection criteria.

    Parameters:
    - cfg (dict): Configuration dictionary.
    - pt_mins (list): List of minimum pt values.
    - pt_maxs (list): List of maximum pt values.
    - bdt_cut_mins (list): List of minimum BDT cut values.
    - bdt_cut_maxs (list): List of maximum BDT cut values.

    Returns:
    - tuple: A tuple containing:
        - dfs_data (list): list of pandas.DataFrame containing the selected data (one per pT bin)
        - dfs_prd_bkg_orig (list): list of lists of pandas.DataFrame containing the selected
            partially reco decays background MC data (one list with all the contributions per pT bin)
        - dfs_prd_bkg_av (list): list of pandas.DataFrame containing the average of selected
            partially reco decays background MC data (one per pT bin)
        - fracs_singlecontr (list): list of lists of fractions for each contribution of correlated bkg
            to be fixed in the fit
    """

    # pt integrated dataframes
    # load data
    df_data = pd.DataFrame()
    for file in cfg["inputs"]["data"]:
        df_data = pd.concat([df_data, pd.read_parquet(file)])
    # load mc
    df_mc = pd.DataFrame()
    for file in cfg["inputs"]["mc"]:
        df_mc = pd.concat([df_mc, pd.read_parquet(file)])

    correlated_bkgs = cfg["correlated_bkgs"]

    dfs_data, dfs_prd_bkg_orig, dfs_prd_bkg_av, fracs_singlecontr = ([] for _ in range(4))
    for ipt, (pt_min, pt_max, bdt_cut_min, bdt_cut_max) in enumerate(zip(pt_mins, pt_maxs,
                                                                         bdt_cut_mins, bdt_cut_maxs)):

        selection_string = f"({pt_min} < fPt < {pt_max} and {bdt_cut_min} < ML_output < {bdt_cut_max})"

        # data
        dfs_data.append(df_data.query(selection_string))

        # mc
        df_mc_pt = df_mc.query(selection_string)
        dfs_mc_sig = df_mc_pt.query("fFlagMcMatchRec == -1 or fFlagMcMatchRec == 1")

        # mc correlated bkgs
        dfs_prd_bkg_orig_pt, fracs_pt = [], []
        den_norm = len(dfs_mc_sig) * cfg["signal_br"]["pdg"] / cfg["signal_br"]["sim"]
        for bkg in correlated_bkgs:
            dfs_prd_bkg_pt = df_mc_pt.query("fFlagMcMatchRec == 8 "
                                            f"and fPdgCodeBeautyMother == {bkg['beauty_id']} and "
                                            f"fPdgCodeCharmMother == {bkg['charm_id']}")
            dfs_prd_bkg_orig_pt.append(dfs_prd_bkg_pt)
            # store all fractions to fix it later in the fit
            fracs_pt.append(len(dfs_prd_bkg_pt) * bkg["br_pdg"] / bkg["br_sim"] / den_norm)
        fracs_singlecontr.append(fracs_pt)
        # sum of all fractions to normalise the weighted average
        sum_fracs_pt = sum(fracs_pt)
        fracs_pt_norm = [frac / sum_fracs_pt for frac in fracs_pt]

        # dfs for each source separately
        dfs_prd_bkg_orig.append(dfs_prd_bkg_orig_pt)

        # dfs for weighted average
        dfs_prd_bkg_sampled = []
        for frac, df_bkg in zip(fracs_pt_norm, dfs_prd_bkg_orig_pt):
            dfs_prd_bkg_sampled.append(df_bkg.sample(frac=frac, random_state=42))

        dfs_prd_bkg_av.append(pd.concat(dfs_prd_bkg_sampled))

    return dfs_data, dfs_prd_bkg_orig, dfs_prd_bkg_av, fracs_singlecontr


def get_fixed_parameter(par, unc, string):
    """
    Return the value of a fixed parameter based on the given string.

    Parameters:
    - par (float): The value of the parameter.
    - unc (float): The uncertainty of the parameter.
    - string (str): The string indicating the type of fixed parameter.

    Returns:
    - float: The value of the fixed parameter.

    Raises:
    - FATAL: If the string does not match any of the patterns.

    """
    if string == "fixed":
        return par
    if string == "fixed_plus_unc":
        return par + unc
    if string == "fixed_minus_unc":
        return par - unc

    # Check if the string matches the pattern fixed_plus_XX_perc
    pattern = r'^fixed_plus_\d+_perc$'
    if re.match(pattern, string):
        # Extract the number
        number = int(string.split('_')[-2])
        return par + par * number / 100

    # Check if the string matches the pattern fixed_minus_XX_perc
    pattern = r'^fixed_minus_\d+_perc$'
    if re.match(pattern, string):
        # Extract the number
        number = int(string.split('_')[-2])
        return par - par * number / 100

    Logger(
        f"The string {string} does not match any of the available patterns", "FATAL"
    )
    return None


def build_data_handlers(trial, df_pt, df_mc_prd_bkg_pt, df_mc_prd_bkg_av_pt):
    """
    Build data handlers for the given trial.

    Parameters:
    - trial (dict): A dictionary containing the trial parameters (mins, maxs, sgn_funcs, bkg_funcs,
        sigma, mean, use_bkg_templ, bincounting_nsigma).
    - df_pt (pandas.DataFrame): The data frame for the data.
    - df_mc_prd_bkg_pt (list): List of pandas.DataFrame for the MC partly reco decays (each contribution separately).
    - df_mc_prd_bkg_av_pt (pandas.DataFrame): The data frame for the MC partly reco decays (average of all contributions).

    Returns:
    - data_hdl (flarefly.DataHandler): The data handler for the data.
    - data_hdl_bkg (list): The list of flarefly.DataHandler for the MC partly reco decays,
        or None if use_bkg_templ is False.
    """
    # data
    data_hdl = DataHandler(df_pt, var_name="fM",
                           limits=[trial["mins"], trial["maxs"]],
                           nbins=round((trial["maxs"]-trial["mins"])/0.01))

    # mc partly reco decays
    data_hdl_bkg = []
    if trial["use_bkg_templ"]:
        if trial["bkg_templ_opt"] == 0:
            for df_mc in df_mc_prd_bkg_pt:
                data_hdl_bkg.append(
                    DataHandler(df_mc, var_name="fM", limits=[trial["mins"], trial["maxs"]],
                                nbins=round((trial["maxs"]-trial["mins"])/0.01)))
        else:
            data_hdl_bkg.append(
                DataHandler(df_mc_prd_bkg_av_pt, var_name="fM", limits=[trial["mins"], trial["maxs"]],
                            nbins=round((trial["maxs"]-trial["mins"])/0.01)))
    else:
        data_hdl_bkg = None

    return data_hdl, data_hdl_bkg


def build_fitter(
        trial, data_hdl, data_hdl_prd_bkg,
        mean_with_unc, sigma_with_unc, fitter_suffix, cfg, fracs
    ):  # pylint: disable=too-many-arguments,too-many-branches,too-many-statements # noqa: 121, 125
    """
    Build a flarefly mass fitter for fitting the data candidate distribution.

    Parameters:
    - trial (dict): A dictionary containing the trial parameters (mins, maxs, sgn_funcs, bkg_funcs,
        bkg_funcs, sigma, mean, use_bkg_templ, bincounting_nsigma).
    - data_hdl (flarefly.DataHandler): The data handler for the data.
    - data_hdl_prd_bkg (list): The list of flarefly.DataHandler for the MC partly reco decays.
    - mean_with_unc (tuple): the mean value of the data with uncertainty.
    - sigma_with_unc (tuple): the sigma value of the data with uncertainty.
    - fitter_suffix (str): A suffix to be added to the fitter name.
    - cfg (dict): config for labels of correlated backgrounds
    - fracs (list): fractions of correlated backgrounds
    - sigma_mc ()

    Returns:
    - fitter (flarefly.F2MassFitter): The mass fitter object.
    """
    correlated_bkgs = cfg["correlated_bkgs"]
    bkg_funcs = trial["bkg_funcs"].copy()
    if not isinstance(bkg_funcs, list):
        bkg_funcs = [bkg_funcs]
    sgn_funcs = trial["sgn_funcs"]
    if not isinstance(sgn_funcs, list):
        sgn_funcs = [sgn_funcs]
    label_bkg_pdf = ["Comb. bkg"]
    if trial["use_bkg_templ"]:
        if trial["bkg_templ_opt"] == 0:
            for i_bkg, bkg in enumerate(correlated_bkgs):
                label_bkg_pdf.insert(i_bkg, bkg["name"])
                bkg_funcs.insert(i_bkg, "kde_grid")
        elif trial["bkg_templ_opt"] == 1:
            label_bkg_pdf.insert(0, "Correlated backgrounds")
            bkg_funcs.insert(0, "kde_grid")

    fitter = F2MassFitter(
        data_hdl,
        sgn_funcs,
        bkg_funcs,
        name=f"b0_{fitter_suffix}",
        label_signal_pdf=[r"$\mathrm{B}^{0}$ signal"],
        label_bkg_pdf=label_bkg_pdf,
        verbosity=0
    )
    if trial["use_bkg_templ"]:
        if trial["bkg_templ_opt"] == 0:
            for i_bkg, (bkg, hdl_bkg) in enumerate(zip(correlated_bkgs, data_hdl_prd_bkg)):
                fitter.set_background_kde(i_bkg, hdl_bkg)
                if i_bkg == 0:
                    if trial["fix_correlated_bkg_to_signal"]:
                        fitter.fix_bkg_frac_to_signal_pdf(i_bkg, 0, fracs[i_bkg])
                    else:
                        fitter.set_background_initpar(i_bkg, "frac", 0.01, limits=[0., 1.])
                else:
                    denom = (
                        data_hdl_prd_bkg[0].get_norm() * correlated_bkgs[0]["br_pdg"] / correlated_bkgs[0]["br_sim"])
                    fitter.fix_bkg_frac_to_bkg_pdf(
                        i_bkg, 0, hdl_bkg.get_norm() * bkg["br_pdg"] / bkg["br_sim"] / denom)

        elif trial["bkg_templ_opt"] == 1:
            fitter.set_background_kde(0, data_hdl_prd_bkg[0])
            if trial["fix_correlated_bkg_to_signal"]:
                fitter.fix_bkg_frac_to_signal_pdf(0, 0, sum(fracs))
            else:
                fitter.set_background_initpar(0, "frac", 0.5, limits=[0., 1.])

    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])

    if "fixed" in trial["mean"]:  # if mean is fixed, set the particle mass
        mean_to_fix = get_fixed_parameter(mean_with_unc[0], mean_with_unc[1], trial["mean"])
        fitter.set_particle_mass(0, mass=mean_to_fix, fix=True)
    else:  # mean is free
        fitter.set_particle_mass(0, pdg_id=511, fix=False)

    if "fixed" in trial["sigma"]:  # if sigma is fixed, set the sigma
        sigma_to_fix = get_fixed_parameter(sigma_with_unc[0], sigma_with_unc[1], trial["sigma"])
        fitter.set_signal_initpar(0, "sigma", sigma_to_fix, fix=True)
    else:  # sigma is free
        fitter.set_signal_initpar(0, "sigma", sigma_with_unc[0], fix=False,
                                  limits=[sigma_with_unc[0] * 0.5, sigma_with_unc[0] * 1.5])

    icombbkg = len(data_hdl_prd_bkg)
    fitter.set_background_initpar(icombbkg, "lam", -1.2, limits=[-10., 10.])
    fitter.set_background_initpar(icombbkg, "c1", -0.05, limits=[-2., 2.])
    fitter.set_background_initpar(icombbkg, "c2", 0.008, limits=[0.000, 0.2])

    return fitter


def fit(fitter, cfg, i_trial, suffix):  # pylint: disable=too-many-locals
    """
    Perform a fit using the given fitter object and configuration.

    Parameters:
    - fitter (flarefly.F2MassFitter): The mass fitter object.
    - cfg (dict): Configuration dictionary.
    - i_trial (int): The trial number.
    - suffix (str): The suffix for the output file.

    Returns:
    - output_dict (dict): A dictionary containing the fit results.
    """
    result = fitter.mass_zfit()
    if result.converged:
        rawy, rawy_unc = fitter.get_raw_yield(0)
        if cfg["multitrial"]["bincounting_nsigma"]:  # if there is at least one nsigma
            rawy_bincounting, rawy_bincounting_unc = zip(
                *[fitter.get_raw_yield_bincounting(0, nsigma=nsigma)
                    for nsigma in cfg["multitrial"]["bincounting_nsigma"]]
            )
        else:
            rawy_bincounting, rawy_bincounting_unc = None, None
        significance, significance_unc = fitter.get_significance(0)
        soverb, soverb_unc = fitter.get_signal_over_background(0)
        mean, mean_unc = fitter.get_mass(0)
        sigma, sigma_unc = fitter.get_sigma(0)
        chi2_ndf = fitter.get_chi2_ndf()
        if cfg["save_all_fits"]:
            fig, _ = fitter.plot_mass_fit(
                style="ATLAS",
                figsize=(8, 8),
                axis_title=r"$M(\mathrm{D^-\pi^+})$ (GeV/$c^2$)"
            )
            if not os.path.exists(os.path.join(cfg["output_dir"], cfg["output_dir_fits"])):
                os.makedirs(os.path.join(cfg["output_dir"], cfg["output_dir_fits"]))
            fig.savefig(
                os.path.join(cfg["output_dir"], cfg["output_dir_fits"], f"mass_fit_{suffix}.pdf")
            )
    else:
        rawy, rawy_unc = None, None
        rawy_bincounting, rawy_bincounting_unc = [None] * len(cfg["multitrial"]["bincounting_nsigma"]), [None] * len(cfg["multitrial"]["bincounting_nsigma"])  # pylint: disable=line-too-long # noqa: 501
        significance, significance_unc = None, None
        soverb, soverb_unc = None, None
        mean, mean_unc = None, None
        sigma, sigma_unc = None, None
        chi2_ndf = None
        Logger(
            f"The fit for the trial {i_trial} did not converge, skipping this trial",
            "WARNING"
        )

    output_dict = {
        "rawy": rawy, "rawy_unc": rawy_unc, "significance": significance,
        "significance_unc": significance_unc, "soverb": soverb, "soverb_unc": soverb_unc,
        "mean": mean, "mean_unc": mean_unc, "sigma": sigma, "sigma_unc": sigma_unc,
        "chi2_ndf": chi2_ndf
    }
    for i_nsigma, nsigma in enumerate(cfg["multitrial"]["bincounting_nsigma"]):
        output_dict[f"rawy_bincounting_{nsigma}"] = rawy_bincounting[i_nsigma]
        output_dict[f"rawy_bincounting_{nsigma}_unc"] = rawy_bincounting_unc[i_nsigma]

    return output_dict


def get_rms_shift_sum_quadrature(df, h_rawy, i_pt, rel=False):
    """
    Calculate the sum in quadrature of the RMS and shift from the central value for raw yields.

    Parameters:
        df (pandas.DataFrame): DataFrame containing the raw yields.
        h_rawy (hist): histogram with central raw yields
        i_pt (int): pt bin index
        rel (bool): If True, return the relative uncertainty.

    Returns:
        float: The sum in quadrature of the RMS and shift from the central value for raw yields.
    """

    central_rawy = h_rawy.values()[i_pt]

    if rel:
        return np.sqrt(
            np.std(df["rawy"])**2 +
            (np.mean(df["rawy"]) - central_rawy)**2
        ) / central_rawy

    return np.sqrt(
        np.std(df["rawy"])**2 +
        (np.mean(df["rawy"]) - central_rawy)**2
    )


def dump_results_to_root(dfs, cfg, h_rawy, cut_set):
    """
    Dump the results to a ROOT file.

    Parameters:
        dfs (list of pandas.DataFrame): List of dataframes containing the data for each pt bin.
        cfg (dict): Configuration dictionary.
        h_rawy (hist): Histogram with central raw yields
        cut_set (dict): Dictionary containing the cut sets.

    Returns:
        None
    """
    pt_mins = cut_set["pt"]["mins"]
    pt_maxs = cut_set["pt"]["maxs"]
    pt_edges = np.asarray(pt_mins + [pt_maxs[-1]], "d")
    pt_bins = cfg["multitrial"]["pt_bins"]
    if pt_bins is None:
        pt_bins = list(range(len(pt_mins)))

    rms_shifts = []
    assigned_syst = []

    cols_to_save = [
        "rawy", "rawy_unc", "significance", "significance_unc",
        "soverb", "soverb_unc", "mean", "mean_unc", "sigma", "sigma_unc", "chi2_ndf"
    ]

    idx_assigned_syst = 0
    for i_pt in range(len(pt_mins)):
        if i_pt not in pt_bins:
            rms_shifts.append(0)
            assigned_syst.append(0)
            continue
        rms_shifts.append(get_rms_shift_sum_quadrature(dfs[i_pt], h_rawy, i_pt, rel=True))
        assigned_syst.append(cfg["assigned_syst"][idx_assigned_syst])
        idx_assigned_syst += 1

    if not os.path.exists(cfg["output_dir"]):
        os.makedirs(cfg["output_dir"])
    output_file_name = os.path.join(cfg["output_dir"], "raw_yields_systematic.root")
    with uproot.recreate(output_file_name) as f:
        for pt_bin in pt_bins:
            suffix = f"_{pt_mins[pt_bin] * 10:.0f}_{pt_maxs[pt_bin] * 10:.0f}"
            f[f"df{suffix}"] = dfs[pt_bin][cols_to_save]

        f["rms_shifts_sum_quadrature"] = (np.array(rms_shifts), pt_edges)
        f["assigned_syst"] = (np.array(assigned_syst), pt_edges)


def process_trial(args):
    """Process a single trial."""
    trial, dfs_data, dfs_mc_prd_bkg, dfs_mc_prd_bkg_av, h_mean_mc, h_sigma_mc, fracs_prd_bkg, cfg, i_pt, pt_min, pt_max, i_trial = args

    suffix = f"{pt_min*10:.0f}_{pt_max*10:.0f}_{i_trial}"
    mean_with_unc = [h_mean_mc.values()[i_pt], h_mean_mc.errors()[i_pt]]
    sigma_with_unc = [h_sigma_mc.values()[i_pt], h_sigma_mc.errors()[i_pt]]

    data_hdl, data_hdl_prd_bkg = build_data_handlers(
        trial, dfs_data[i_pt], dfs_mc_prd_bkg[i_pt], dfs_mc_prd_bkg_av[i_pt]
    )

    fitter = build_fitter(
        trial, data_hdl, data_hdl_prd_bkg, mean_with_unc, sigma_with_unc,
        suffix, cfg, fracs_prd_bkg[i_pt]
    )

    trial_dict = fit(fitter, cfg, i_trial, suffix)
    trial_renamed = trial.copy()
    trial_renamed['sigma_type'] = trial_renamed.pop('sigma')
    trial_renamed['mean_type'] = trial_renamed.pop('mean')
    trial_dict.update(trial_renamed)

    return trial_dict


def multi_trial(config_file_name: str, draw_only: bool = False):
    """
    Perform multiple trials based on the given configuration file in parallel.

    Parameters:
        config_file_name (str): The path to the configuration file.
    Returns:
        None
    """

    # Load configuration and reference results
    with open(config_file_name, "r", encoding="utf-8") as file_config:
        cfg = yaml.safe_load(file_config)

    zfit.run.set_cpus_explicit(
        intra=cfg["multiprocessing"]["zfit_cpus"]["intra"],
        inter=cfg["multiprocessing"]["zfit_cpus"]["inter"],
    )

    with uproot.open(cfg["reference_fits"]) as file_ref:
        h_rawy = file_ref["h_rawyields"]
        h_sigma = file_ref["h_sigmas"]
        h_mean_mc = file_ref["h_means_mc"]
        h_sigma_mc = file_ref["h_sigmas_mc"]

    # Load cut set
    with open(cfg["inputs"]["cutset"], "r", encoding="utf-8") as file_cutset:
        cut_set = yaml.safe_load(file_cutset)

    pt_mins = cut_set["pt"]["mins"]
    pt_maxs = cut_set["pt"]["maxs"]
    bdt_cut_mins = cut_set["ML_output"]["mins"]
    bdt_cut_maxs = cut_set["ML_output"]["maxs"]

    # Check if only one signal function is used per trial
    multitrial_cfg = cfg["multitrial"]
    for signal_func in multitrial_cfg["sgn_funcs"]:
        if isinstance(signal_func, list) and len(signal_func) > 1:
            Logger("Only one signal function is allowed per trial", "FATAL")

    # Define all the trials
    trials = list(itertools.product(*(multitrial_cfg[var] for var in MULTITRIAL_PARAMS)))
    trials = [dict(zip(MULTITRIAL_PARAMS, trial)) for trial in trials]

    # Get input data
    dfs_data, dfs_mc_prd_bkg, dfs_mc_prd_bkg_av, fracs_prd_bkg = get_input_data(
        cfg, pt_mins, pt_maxs, bdt_cut_mins, bdt_cut_maxs
    )

    dfs = []
    idx_assigned_syst = 0
    print(pt_mins, pt_maxs)

    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        if multitrial_cfg["pt_bins"] is not None and i_pt not in multitrial_cfg["pt_bins"]:
            dfs.append(None)
            continue
        if not draw_only:

            # Prepare arguments for parallel execution of trials
            args = [
                (trial, dfs_data, dfs_mc_prd_bkg, dfs_mc_prd_bkg_av, h_mean_mc, h_sigma_mc, fracs_prd_bkg, cfg, i_pt, pt_min, pt_max, i_trial)
                for i_trial, trial in enumerate(trials)
            ]

            # Parallelize the trials
            with ProcessPoolExecutor(max_workers=cfg["multiprocessing"]["max_workers"]) as executor:
                trial_results = list(executor.map(process_trial, args))

            # Save results
            df_trials = pd.DataFrame(trial_results)
            if not os.path.exists(cfg["output_dir"]):
                os.makedirs(cfg["output_dir"])
            df_trials.to_parquet(
                os.path.join(cfg["output_dir"], f"raw_yields_{pt_min*10:.0f}_{pt_max*10:.0f}.parquet")
            )
        else:
            df_trials = pd.read_parquet(
                os.path.join(cfg["output_dir"], f"raw_yields_{pt_min*10:.0f}_{pt_max*10:.0f}.parquet")
            )
        dfs.append(df_trials)

        # Draw results
        draw_multitrial(df_trials, cfg, pt_min, pt_max, idx_assigned_syst, h_rawy, h_sigma)
        idx_assigned_syst += 1

    dump_results_to_root(dfs, cfg, h_rawy, cut_set)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Multitrial B mesons')
    parser.add_argument('--configFile', '-c', metavar='text',
                        help='Path to the configuration file')
    parser.add_argument('--draw-only', action='store_true',
                        help='Draw only the results')
    args = parser.parse_args()

    # "bincounting_nsigma" removed so that we do not fit multiple times for each nsigma
    MULTITRIAL_PARAMS = ["mins", "maxs", "sgn_funcs", "bkg_funcs", "sigma", "mean",
                         "use_bkg_templ", "bkg_templ_opt", "fix_correlated_bkg_to_signal"]
    multi_trial(args.configFile, args.draw_only)

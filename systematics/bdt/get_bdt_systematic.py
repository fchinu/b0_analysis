"""
Perform evaluation of systematic uncertainties on BDT cuts.

Run: python get_bdt_systematic.py config_file
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
import uproot  # noqa; E402
import yaml  # noqa; E402
import matplotlib.pyplot as plt  # noqa; E402
from matplotlib.patches import Rectangle  # noqa; E402
from PyPDF2 import PdfMerger
import re
import numpy as np  # noqa; E402
import pandas as pd  # noqa; E402
import zfit
from flarefly.utils import Logger  # noqa; E402
from flarefly.data_handler import DataHandler  # noqa; E402
from flarefly.fitter import F2MassFitter  # noqa; E402


def get_axis_range(df, column, central_value, central_unc, is_ratio=False):
    """
    Get the maximum value for the y-axis in a plot.

    Parameters:
        - df (pandas.DataFrame): The dataframe containing the data.
        - column (str): The column name to consider.
        - central_value (float): The central value.
        - central_unc (float): The central uncertainty.
        - is_ratio (bool): Boolean indicating whether the plot is a ratio plot.
    Returns:
        - float: The maximum value for the y-axis.
    """
    if is_ratio:
        max_value = max(
            df[column].max() / central_value,
            1 + central_unc / central_value,
        )
        min_value = min(
            df[column].min() / central_value,
            1 - central_unc / central_value,
        )
    else:
        max_value = max(
            df[column].max(),
            central_value + central_unc,
        )
        min_value = min(
            df[column].min(),
            central_value - central_unc,
        )
    return min_value * 0.9, max_value * 1.1


def draw_cut_variation(df, config, pt_min, pt_max, idx_assigned_syst):  # pylint: disable=too-many-locals, too-many-statements # noqa: E501
    """
    Draws the results of the cut variaiton.

    Parameters:
    - df (pandas.DataFrame): DataFrame containing the data to be plotted.
    - config (dict): Configuration dictionary.
    - pt_min (float): Minimum pT value for the plot title and output filename.
    - pt_max (float): Maximum pT value for the plot title and output filename.
    - idx_assigned_syst (int): Index of the assigned systematic uncertainty.
    """
    # Create figures
    fig, axs = plt.subplots(2, 3, figsize=(32, 20))

    # get the central values
    with uproot.open(config["fit"]["fit_file"], encoding="utf8") as f:
        h_rawy = f["h_rawyields"]
        i_pt = np.digitize((pt_min + pt_max) / 2, h_rawy.axis().edges()) - 1

        central_rawy = f["h_rawyields"].values()[i_pt]
        central_rawy_unc = f["h_rawyields"].errors()[i_pt]
        central_significance = f["h_significance"].values()[i_pt]
        central_significance_unc = f["h_significance"].errors()[i_pt]
        central_s_over_b = f["h_soverb"].values()[i_pt]
        central_s_over_b_unc = f["h_soverb"].errors()[i_pt]

    with uproot.open(config["efficiency_file"], encoding="utf8") as f:
        central_efficiency = f["h_eff"].values()[i_pt]
        central_efficiency_unc = f["h_eff"].errors()[i_pt]
    # Get the corrected raw yields
    central_dict = {
        "rawy": central_rawy,
        "rawy_unc": central_rawy_unc,
        "eff": central_efficiency,
        "eff_unc": central_efficiency_unc,
    }
    central_corr_rawy, central_corr_rawy_unc = get_corr_rawy(central_dict)

    df = df.dropna()  # Drop rows with None values (not converged fits)
    # Raw yields
    axs[0, 0].scatter(range(1, len(df) + 1), df["rawy"] / central_rawy, s=100)
    axs[0, 0].set_title('Raw Yields', fontsize=20)
    axs[0, 0].set_ylabel('Raw Yield/Raw Yield (central)', fontsize=16)
    axs[0, 0].set_xlabel('Cut set', fontsize=16)
    # Draw the central values
    axs[0, 0].axhline(y=1, color='r', linestyle='--')
    axs[0, 0].add_patch(
        Rectangle(
            (0, 1 - central_rawy_unc/central_rawy),
            len(df) + 1, 2 * central_rawy_unc/central_rawy,
            color='r', alpha=0.3, zorder=0,
            label=r'central $\pm \sigma$'
        )
    )
    axs[0, 0].set_ylim(*get_axis_range(df, "rawy", central_rawy, central_rawy_unc, is_ratio=True))
    axs[0, 0].legend(fontsize=12)

    # Efficiencies
    axs[0, 1].scatter(range(1, len(df) + 1), df["eff"] / central_efficiency, s=100)
    axs[0, 1].set_title('Efficiency', fontsize=20)
    axs[0, 1].set_ylabel('Efficiency/Efficiency (central)', fontsize=16)
    axs[0, 1].set_xlabel('Cut set', fontsize=16)
    # Draw the central values
    axs[0, 1].axhline(y=1, color='r', linestyle='--')
    axs[0, 1].add_patch(
        Rectangle(
            (0, 1 - central_efficiency_unc/central_efficiency),
            len(df) + 1, 2 * central_efficiency_unc/central_efficiency,
            color='r', alpha=0.3, zorder=0,
            label=r'central $\pm \sigma$'
        )
    )
    axs[0, 1].set_ylim(*get_axis_range(
        df, "eff", central_efficiency,
        central_efficiency_unc, is_ratio=True
    ))
    axs[0, 1].legend(fontsize=12)

    # Significance
    axs[0, 2].errorbar(
        range(1, len(df) + 1), df["significance"],
        df["significance_unc"], fmt='o', markersize=10
    )
    axs[0, 2].set_title('Significance', fontsize=20)
    axs[0, 2].set_ylabel('Significance', fontsize=16)
    axs[0, 2].set_xlabel('Cut set', fontsize=16)
    # Draw the central values
    axs[0, 2].axhline(y=central_significance, color='r', linestyle='--')
    axs[0, 2].add_patch(
        Rectangle(
            (0, central_significance - central_significance_unc),
            len(df) + 1, 2 * central_significance_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'central $\pm \sigma$'
        )
    )
    axs[0, 2].set_ylim(*get_axis_range(
        df, "significance", central_significance,
        central_significance_unc
    ))
    axs[0, 2].legend(fontsize=12)

    # S/B
    axs[1, 0].errorbar(
        range(1, len(df) + 1), df["soverb"],
        df["soverb_unc"], fmt='o', markersize=10
    )
    axs[1, 0].set_title('S/B', fontsize=20)
    axs[1, 0].set_ylabel('S/B', fontsize=16)
    axs[1, 0].set_xlabel('Cut set', fontsize=16)
    # Draw the central values
    axs[1, 0].axhline(y=central_s_over_b, color='r', linestyle='--')
    axs[1, 0].add_patch(
        Rectangle(
            (0, central_s_over_b - central_s_over_b_unc),
            len(df) + 1, 2 * central_s_over_b_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'central $\pm \sigma$'
        )
    )
    axs[1, 0].set_ylim(*get_axis_range(df, "soverb", central_s_over_b, central_s_over_b_unc))
    axs[1, 0].legend(fontsize=12)

    # Corrected yields
    axs[1, 1].errorbar(
        range(1, len(df) + 1), df["corr_rawy"],
        df["corr_rawy_unc"], fmt='o', markersize=10
    )
    axs[1, 1].set_title(r'$N^\mathrm{raw}/\varepsilon$', fontsize=20)
    axs[1, 1].set_ylabel(r'$N^\mathrm{raw}/\varepsilon$', fontsize=16)
    axs[1, 1].set_xlabel('Cut set', fontsize=16)
    # Draw the central values
    axs[1, 1].axhline(y=central_corr_rawy, color='r', linestyle='--')
    axs[1, 1].add_patch(
        Rectangle(
            (0, central_corr_rawy - central_corr_rawy_unc),
            len(df) + 1, 2 * central_corr_rawy_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'central $\pm \sigma$'
        )
    )
    axs[1, 1].set_ylim(*get_axis_range(df, "corr_rawy", central_corr_rawy, central_corr_rawy_unc))
    axs[1, 1].legend(fontsize=12)

    # Corrected yields histogram
    axs[1, 2].hist(df["corr_rawy"]/central_corr_rawy, bins=20, histtype='step', color='k')
    axs[1, 2].set_title('Ratio', fontsize=20)
    axs[1, 2].set_ylabel('Counts', fontsize=16)
    axs[1, 2].set_xlabel('Ratio/Ratio (central)', fontsize=16)
    # Draw the rms + shift from mean
    rms_shift = get_rms_shift_sum_quadrature(df, config, i_pt, rel=True)
    axs[1, 2].axvline(x=1, color='r', linestyle='--')
    axs[1, 2].add_patch(
        Rectangle(
            (1 - rms_shift, 0),
            2 * rms_shift, np.sqrt(len(df)),
            color='r', alpha=0.3, zorder=0,
            label=r'$\sqrt{\mathrm{RMS}^2 + \Delta^2}$'
        )
    )
    # Draw the assigned systematic uncertainty
    syst_unc = config["assigned_syst"][idx_assigned_syst]
    axs[1, 2].add_patch(
        Rectangle(
            (1 - syst_unc, 0),
            2 * syst_unc, np.sqrt(len(df)),
            color='tab:blue', alpha=0.3, zorder=0,
            label='Assigned syst.'
        )
    )
    x_min = min(
        df["corr_rawy"].min() / central_corr_rawy,
        1 - rms_shift / central_corr_rawy,
        1 - syst_unc
    )
    x_max = max(
        df["corr_rawy"].max() / central_corr_rawy,
        1 + rms_shift / central_corr_rawy,
        1 + syst_unc
    )
    axs[1, 2].set_xlim(x_min * 0.9, x_max * 1.1)
    axs[1, 2].legend(fontsize=12)

    # Save the figure
    if not os.path.exists(os.path.expanduser(f"{config['output']['output_dir']}")):
        os.makedirs(os.path.expanduser(f"{config['output']['output_dir']}"))
    fig.savefig(os.path.join(
        os.path.expanduser(f"{config['output']['output_dir']}"),
        f"BDT_{pt_min * 10:.0f}_{pt_max * 10:.0f}.png"),
        bbox_inches='tight'
    )


def load_data_mc_df(config):
    """
    Load the data from the input file.

    Parameters:
        - config (dict): The configuration object containing input information.
    Returns:
        - df_data (pandas.DataFrame): The data dataframe.
        - df_mc (pandas.DataFrame): The MC dataframe.
    """
    df_data = pd.concat([pd.read_parquet(input_file) for input_file in config["inputs"]["data"]])
    df_mc = pd.concat([pd.read_parquet(input_file) for input_file in config["inputs"]["mc"]])
    if config["inputs"]["mc_for_efficiency"] is not None:
        df_mc_for_efficiency = pd.concat([pd.read_parquet(input_file) for input_file in config["inputs"]["mc_for_efficiency"]])
    else:
        df_mc_for_efficiency = df_mc
    return df_data, df_mc, df_mc_for_efficiency


def get_cuts(config, central_cutset):  # pylint: disable=too-many-locals
    """
    Generate systematic variations of cuts based on the provided configuration and central cutset.

    Parameters:
        - config (dict): Configuration dictionary containing cut variation parameters.
        - central_cutset (dict): Dictionary containing central cut values.

    Returns:
        - min_selections (list of list of float): List of lists containing minimum
            ML output variations for each pt bin.
        - max_selections (list of list of float): List of lists containing maximum
            ML output variations for each pt bin.
    """
    ml_mins = central_cutset["ML_output"]["mins"]
    ml_maxs = central_cutset["ML_output"]["maxs"]

    min_selections, max_selections = [], []
    for i_pt, (ml_min, ml_max) in enumerate(zip(ml_mins, ml_maxs)):
        n_cuts_pos = config["cut_variations"]["n_cuts_pos"][i_pt]
        n_cuts_neg = config["cut_variations"]["n_cuts_neg"][i_pt]
        ml_min_var = config["cut_variations"]["mins"][i_pt]
        ml_max_var = config["cut_variations"]["maxs"][i_pt]
        if config["cut_variations"]["edge"] == "min":
            min_selections.append(
                np.linspace(ml_min_var, ml_mins[i_pt], n_cuts_neg).tolist() +
                np.linspace(ml_mins[i_pt], ml_max_var, n_cuts_pos).tolist()[1:]
            )
            max_selections.append([ml_max] * len(min_selections[-1]))
        elif config["cut_variations"]["edge"] == "max":
            max_selections.append(np.linspace(ml_min_var, ml_max_var, n_cuts).tolist())
            min_selections.append([ml_min] * len(max_selections[-1]))
        else:
            Logger("Invalid edge type.", "FATAL")

    return min_selections, max_selections


def get_fit_config(config, i_pt):
    """
    Extract and construct the fit configuration for a given pt index.

    Parameters:
        - config (dict): Configuration dictionary.
        - i_pt (int): Index for the pt bin.
    Returns:
        - dict_fit_config (dict): A dictionary containing the fit configuration
            with the following keys:
            - "signal_funcs": Signal functions for the given pt bin.
            - "bkg_funcs": Background functions for the given pt bin.
            - "mass_limits": Mass limits for the given pt bin.
            - "use_bkg_templ": Boolean indicating whether to use background templates
                for the given pt bin.
            - "n_bins": Number of bins for the given pt bin.
            - "sigma": Fixed sigma value if specified in the config, otherwise None.
            - "mean": Fixed mean value if specified in the config, otherwise None.
    """
    with open(config["fit"]["fit_config"], 'r', encoding="utf8") as f:
        fit_config = yaml.safe_load(f)

    sigma, mean = None, None
    if config["fit"]["fix_sigma"] or config["fit"]["fix_mean"]:
        with uproot.open(config["fit"]["fit_file"], encoding="utf8") as f:
            if config["fit"]["fix_sigma"]:
                sigma = f["h_sigmas"].values()[i_pt]
            if config["fit"]["fix_mean"]:
                mean = f["h_means"].values()[i_pt]

    dict_fit_config = {
        "signal_funcs": fit_config["fit_configs"]["signal_funcs"][i_pt],
        "bkg_funcs": fit_config["fit_configs"]["bkg_funcs"][i_pt],
        "mass_limits": fit_config["fit_configs"]["mass_limits"][i_pt],
        "use_bkg_templ": fit_config["fit_configs"]["use_bkg_templ"][i_pt],
        "n_bins": fit_config["plot_style"]["n_bins"][i_pt],
        "sigma": sigma,
        "mean": mean,
        "correlated_bkgs": fit_config["fit_configs"]["correlated_bkgs"],
        "signal_br": fit_config["fit_configs"]["signal_br"],
        "corr_bkg_frac": None,
        "i_pt": i_pt
    }

    return dict_fit_config


def build_fitter(df_data, df_mc_prd_bkg, fit_config):
    """
    Build and configure a mass fitter for B0 analysis.

    Parameters:
    - df_data (pd.DataFrame): DataFrame containing the data sample.
    - df_mc_prd_bkg (pd.DataFrame): DataFrame containing the MC
        partly reconstructed background sample.
    - fit_config (dict): Configuration dictionary for the fit.

    Returns:
    - fitter (flarefly.F2MassFitter): Mass fitter.
    """
    # data
    data_hdl = DataHandler(
        df_data, var_name="fM",
        limits=[fit_config["mass_limits"][0], fit_config["mass_limits"][1]],
        nbins=round((fit_config["mass_limits"][1]-fit_config["mass_limits"][0])/0.01)
    )

    # mc partly reco decays
    data_hdl_bkg = DataHandler(
        df_mc_prd_bkg, var_name="fM",
        limits=[fit_config["mass_limits"][0], fit_config["mass_limits"][1]],
        nbins=round((fit_config["mass_limits"][1]-fit_config["mass_limits"][0])/0.01)
    )

    bkg_funcs = fit_config["bkg_funcs"]
    label_bkg_pdf = ["Comb. bkg"]
    if fit_config["use_bkg_templ"]:
        bkg_funcs.insert(0, "kde_grid")
        label_bkg_pdf.insert(0, "Partly reco decays")

    fitter = F2MassFitter(
        data_hdl,
        fit_config["signal_funcs"],
        bkg_funcs,
        name=f"b0_{fit_config['i_pt']}_{fit_config['i_var']}",
        label_signal_pdf=[r"$\mathrm{B}^{0}$ signal"],
        label_bkg_pdf=label_bkg_pdf,
        verbosity=0
    )

    if fit_config["use_bkg_templ"]:
        fitter.set_background_kde(0, data_hdl_bkg)
        fitter.set_background_initpar(0, "frac", 0.01, limits=[0., 1.])

    if fit_config["mean"] is not None:
        fitter.set_signal_initpar(0, "mean", fit_config["mean"], fix=True)
    else:
        fitter.set_particle_mass(0, pdg_id=511, fix=False)
    if fit_config["sigma"] is not None:
        fitter.set_signal_initpar(0, "sigma", fit_config["sigma"], fix=True)
    else:
        fitter.set_signal_initpar(0, "sigma", 0.04, limits=[0.01, 0.08])
    fitter.set_signal_initpar(0, "frac", 0.05, limits=[0., 1.])

    fitter.set_background_initpar(1, "c0", 1.)
    fitter.set_background_initpar(1, "c1", -0.05, limits=[-2, 2.])
    fitter.set_background_initpar(1, "c2", 0.008, limits=[0.000, 0.5])

    fitter.fix_bkg_frac_to_signal_pdf(0, 0, fit_config["corr_bkg_frac"])

    return fitter


def get_fit_results(fitter, result):
    """
    Extract and return the fit results from a fitter object.

    Parameters:
    fitter (flarefly.F2MassFitter): The fitter.
    result (zfit.minimizers.fitresult.FitResult): The fit result.

    Returns:
    dict: A dictionary containing the following keys:
        - "rawy": The raw yield from the fitter (list of two elements if not converged).
        - "significance": The significance of the fit (list of two elements if not converged).
        - "soverb": The signal over background ratio (list of two elements if not converged).
        - "mean": The mean of the signal parameter (None if not converged).
        - "sigma": The sigma of the signal parameter (None if not converged).
        - "chi2": The chi-squared value of the fit (None if not converged).
    """
    if result.converged:
        rawy, rawy_unc = fitter.get_raw_yield(0)
        significance, significance_unc = fitter.get_significance(0)
        soverb, soverb_unc = fitter.get_signal_over_background(0)
        mean, mean_unc = fitter.get_signal_parameter(0, "mu")
        sigma, sigma_unc = fitter.get_signal_parameter(0, "sigma")

        result = {
            "rawy": rawy,
            "rawy_unc": rawy_unc,
            "significance": significance,
            "significance_unc": significance_unc,
            "soverb": soverb,
            "soverb_unc": soverb_unc,
            "mean": mean,
            "mean_unc": mean_unc,
            "sigma": sigma,
            "sigma_unc": sigma_unc,
            "chi2": float(fitter.get_chi2())
        }
    else:
        result = {
            "rawy": None,
            "rawy_unc": None,
            "significance": None,
            "significance_unc": None,
            "soverb": None,
            "soverb_unc": None,
            "mean": None,
            "mean_unc": None,
            "sigma": None,
            "sigma_unc": None,
            "chi2": None
        }
    return result


def get_efficiency(df_mc, df_mc_sel, config, fit_config):
    """
    Calculate the efficiency and its uncertainty.

    Parameters:
    df_mc (pandas.DataFrame): DataFrame containing the MC data.
    df_mc_sel (pandas.DataFrame): DataFrame containing the selected MC data.
    config (dict): Cut variation onfiguration dictionary.
    fit_config (dict): Configuration dictionary containing the fit index.

    Returns:
    tuple: A tuple containing:
        - eff (float): The calculated efficiency.
        - eff_unc (float): The uncertainty of the calculated efficiency.
    """
    with uproot.open(config["efficiency_file"], encoding="utf8") as f:
        trigger_eff = f["h_eff_trigger"].values()[fit_config["i_pt"]]
    eff = len(df_mc_sel) / len(df_mc) * trigger_eff
    # assume trigger efficiency uncertainty is negligible
    eff_unc = np.sqrt(eff * (1 - eff) / len(df_mc))
    return eff, eff_unc


def get_corr_rawy(results):
    """
    Calculate the corrected raw yield and its uncertainty.

    Parameters:
    results (dict): A dictionary containing the following keys:
        - "rawy": The raw yield.
        - "eff": The efficiency.
        - "rawy_unc": The uncertainty in the raw yield.
        - "eff_unc": The uncertainty in the efficiency.

    Returns:
    tuple: A tuple containing:
        - corr_rawy (float): The corrected raw yield.
        - corr_rawy_unc (float): The uncertainty in the corrected raw yield.
    """
    corr_rawy = results["rawy"] / results["eff"]
    corr_rawy_unc = np.sqrt(
        (results["rawy_unc"] / results["rawy"])**2 + (results["eff_unc"] / results["eff"])**2
    ) * corr_rawy
    return corr_rawy, corr_rawy_unc


def run_variation(df_data, df_mc, df_mc_eff, selection, config, fit_config):  # pylint: disable=too-many-locals
    """
    Run the variation for the given selection.

    Args:
        - df_data (pandas.DataFrame): The data dataframe.
        - df_mc (pandas.DataFrame): The MC dataframe.
        - selection (str): The selection string to apply.
    Returns:
        - results (dict): A dictionary containing the results of the variation.
    """
    df_data_sel = df_data.query(selection)
    df_mc_sig = df_mc.query("abs(fFlagMcMatchRec) == 1")
    df_mc_sel = df_mc.query(selection)
    df_mc_sel_sig = df_mc_sel.query("abs(fFlagMcMatchRec) == 1")
    df_mc_eff_sig = df_mc_eff.query("abs(fFlagMcMatchRec) == 1")
    df_mc_eff_sel_sig = df_mc_eff_sig.query(selection)
    df_mc_prd_bkg = df_mc_sel.query("fFlagMcMatchRec == 8")

    dfs_prd_bkg_orig = []
    fracs = []
    den_norm = len(df_mc_sel_sig) * fit_config["signal_br"]["pdg"] / fit_config["signal_br"]["sim"]
    for bkg in fit_config["correlated_bkgs"]:
        df_prd_bkg = df_mc_prd_bkg.query(f"fPdgCodeBeautyMother == {bkg['beauty_id']} and "
                                         f"fPdgCodeCharmMother == {bkg['charm_id']}")
        dfs_prd_bkg_orig.append(df_prd_bkg)
        fracs.append(len(df_prd_bkg) * bkg["br_pdg"] / bkg["br_sim"] / den_norm)
    sum_fracs = sum(fracs)
    fracs_norm = [frac / sum_fracs for frac in fracs]

    dfs_prd_bkg_sampled = []
    for frac, df_bkg in zip(fracs_norm, dfs_prd_bkg_orig):
        dfs_prd_bkg_sampled.append(df_bkg.sample(frac=frac, random_state=42))
    df_prd_bkg_sampled = pd.concat(dfs_prd_bkg_sampled)
    fit_config.update({"corr_bkg_frac": sum(fracs)})

    # get the raw yields
    fitter = build_fitter(df_data_sel, df_prd_bkg_sampled, fit_config)

    results = fitter.mass_zfit()
    if results.converged and config["output"]["save_all_fits"]:
        suffix = f"_{fit_config['i_pt']}_{fit_config['i_var']}"
        out_dir = os.path.join(
            os.path.expanduser(f"{config['output']['output_dir']}"),
            config["output"]["output_dir_fits"]
        )
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        fig, _ = fitter.plot_mass_fit(
            style="ATLAS",
            show_extra_info=True,
            figsize=(8, 8),
            axis_title=r"$M(\mathrm{D^-\pi^+})$ (GeV/$c^2$)"
        )
        fig.savefig(
            os.path.join(out_dir, f"mass_fit_{suffix}.pdf")
        )

    variation_results = get_fit_results(fitter, results)
    eff, eff_unc = get_efficiency(df_mc_eff_sig, df_mc_eff_sel_sig, config, fit_config)
    variation_results.update({"eff": eff, "eff_unc": eff_unc})
    corr_rawy, corr_rawy_unc = get_corr_rawy(variation_results)
    variation_results.update({"corr_rawy": corr_rawy, "corr_rawy_unc": corr_rawy_unc})
    variation_results.update({
        "min_selection": fit_config["min_selection"],
        "max_selection": fit_config["max_selection"]
    })
    del df_data_sel, df_mc_sel, df_mc_sel_sig, df_mc_prd_bkg, fitter
    return variation_results


def get_rms_shift_sum_quadrature(df, cfg, i_pt, rel=False):
    """
    Calculate the sum in quadrature of the RMS and shift from the central value for raw yields.

    Parameters:
        df (pandas.DataFrame): DataFrame containing the raw yields.
        cfg (dict): Configuration dictionary.
        i_pt (int): Index of pt.
        rel (bool): If True, return the relative uncertainty.

    Returns:
        float: The sum in quadrature of the RMS and shift from the central value for raw yields.
    """
    with uproot.open(cfg["fit"]["fit_file"]) as f:
        h_rawy = f["h_rawyields"]
    central_rawy = h_rawy.values()[i_pt]
    with uproot.open(cfg["efficiency_file"]) as f:
        h_eff = f["h_eff"]
    central_eff = h_eff.values()[i_pt]
    central_corry = central_rawy / central_eff

    if rel:
        return np.sqrt(
            np.std(df["corr_rawy"])**2 +
            (np.mean(df["corr_rawy"]) - central_corry)**2
        ) / central_corry

    return np.sqrt(
        np.std(df["corr_rawy"])**2 +
        (np.mean(df["corr_rawy"]) - central_corry)**2
    )


def dump_results_to_root(dfs, cfg, cut_set):
    """
    Dump the results to a ROOT file.

    Parameters:
        dfs (list of pandas.DataFrame): List of dataframes containing the data for each pt bin.
        cfg (dict): Configuration dictionary.
        cut_set (dict): Dictionary containing the cut sets.

    Returns:corr_rawy
        None
    """
    pt_mins = cut_set["pt"]["mins"]
    pt_maxs = cut_set["pt"]["maxs"]
    pt_edges = np.asarray(pt_mins + [pt_maxs[-1]], "d")
    pt_bins = cfg["cut_variations"]["pt_bins"]
    if pt_bins is None:
        pt_bins = list(range(len(pt_mins)))

    rms_shifts = []
    assigned_syst = []

    idx_assigned_syst = 0
    for i_pt in range(len(pt_mins)):
        if i_pt not in pt_bins:
            rms_shifts.append(0)
            assigned_syst.append(0)
            continue
        rms_shifts.append(get_rms_shift_sum_quadrature(dfs[i_pt], cfg, i_pt, rel=True))
        assigned_syst.append(cfg["assigned_syst"][idx_assigned_syst])
        idx_assigned_syst += 1

    if not os.path.exists(cfg["output"]["output_dir"]):
        os.makedirs(cfg["output"]["output_dir"])
    output_file_name = os.path.join(cfg["output"]["output_dir"], "bdt_systematic.root")
    with uproot.recreate(output_file_name) as f:
        for pt_bin in pt_bins:
            suffix = f"_{pt_mins[pt_bin] * 10:.0f}_{pt_maxs[pt_bin] * 10:.0f}"
            f[f"df{suffix}"] = dfs[pt_bin]

        f["rms_shifts_sum_quadrature"] = (np.array(rms_shifts), pt_edges)
        f["assigned_syst"] = (np.array(assigned_syst), pt_edges)


def merge_and_clean_pdfs(config, pt_min, pt_max, i_pt):
    """
    Merge all PDFs for a given pt bin into a single file and delete individual files.

    Args:
        config (dict): Configuration dictionary.
        pt_min (float): Minimum pt value for the bin.
        pt_max (float): Maximum pt value for the bin.
        i_pt (int): Index of the pt bin.
    """
    pt_suffix = f"{pt_min * 10:.0f}_{pt_max * 10:.0f}"
    pdf_dir = os.path.join(
            os.path.expanduser(f"{config['output']['output_dir']}"),
            config["output"]["output_dir_fits"]
        )
    merged_pdf_path = os.path.join(pdf_dir, f"mass_fits_{pt_suffix}_merged.pdf")
    
    pdf_merger = PdfMerger()

    def extract_sort_key(filename):
        """
        Extract the numerical parts from the filename as a tuple for sorting.

        Args:
            filename (str): The filename to process.

        Returns:
            tuple: A tuple of numbers extracted from the filename.
        """
        match = re.findall(r'(\d+)', filename)
        return tuple(map(int, match)) if match else (0,)
    
    pdf_files = sorted([os.path.join(pdf_dir, f) for f in os.listdir(pdf_dir) if f.startswith(f"mass_fit__{i_pt}")], key=extract_sort_key)

    for pdf in pdf_files:
       pdf_merger.append(pdf)

    pdf_merger.write(merged_pdf_path)
    pdf_merger.close()

    # Delete individual files
    for pdf in pdf_files:
       os.remove(pdf)

def cut_variation(config_file_name, draw_only=False):  # pylint: disable=too-many-locals
    """
    Perform systematic variations on BDT cuts and save the results.

    Parameters:
        config_file_name (str): Path to the configuration file in YAML format.
    """
    with open(config_file_name, 'r', encoding="utf8") as f:
        config = yaml.safe_load(f)

    with open(config["central_cutset"], 'r', encoding="utf8") as f:
        central_cutset = yaml.safe_load(f)

    pt_mins = central_cutset["pt"]["mins"]
    pt_maxs = central_cutset["pt"]["maxs"]

    min_selections, max_selections = get_cuts(config, central_cutset)
    print(min_selections)
    print(max_selections)

    df_data, df_mc, df_mc_eff = load_data_mc_df(config)

    idx_assigned_syst = 0
    out_dfs = []
    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        if not draw_only:
            if config["cut_variations"]["pt_bins"] is not None:
                if i_pt not in config["cut_variations"]["pt_bins"]:
                    out_dfs.append(None)
                    continue
            fit_config = get_fit_config(config, i_pt)
            df_data_pt = df_data.query(f"{pt_min} < fPt < {pt_max}")
            df_mc_pt = df_mc.query(f"{pt_min} < fPt < {pt_max}")
            df_mc_eff_pt = df_mc_eff.query(f"{pt_min} < fPt < {pt_max}")
            results = []
            with ProcessPoolExecutor(max_workers=config["max_workers"]) as executor:
                for i_var, (min_selection, max_selection) in enumerate(zip(min_selections[i_pt], max_selections[i_pt])):  # pylint: disable=line-too-long # noqa: E501
                    fit_config.update({
                        "i_var": i_var,
                        "min_selection": min_selection,
                        "max_selection": max_selection
                    })
                    selection = f"{min_selection} < ML_output < {max_selection}"
                    results.append(executor.submit(
                        run_variation, df_data_pt, df_mc_pt, df_mc_eff_pt,
                        selection, config, fit_config.copy()
                    ))

            merge_and_clean_pdfs(config, pt_min, pt_max, i_pt)
            out_df = []
            for result in results:
                result = result.result()
                # Wrap into list to avoid ValueError: If using all scalar values, you must pass an index
                out_df.append(pd.DataFrame([result]))

            out_df = pd.concat(out_df)
            if not os.path.exists(os.path.expanduser(f"{config['output']['output_dir']}")):
                os.makedirs(os.path.expanduser(f"{config['output']['output_dir']}"))
            out_df.to_parquet(os.path.join(
                os.path.expanduser(f"{config['output']['output_dir']}"),
                f"BDT_{pt_min * 10:.0f}_{pt_max * 10:.0f}.parquet")
            )
        else:
            out_df = pd.read_parquet(os.path.join(
                os.path.expanduser(f"{config['output']['output_dir']}"),
                f"BDT_{pt_min * 10:.0f}_{pt_max * 10:.0f}.parquet")
            )
        out_dfs.append(out_df)

        draw_cut_variation(out_df, config, pt_min, pt_max, idx_assigned_syst)
        idx_assigned_syst += 1
    dump_results_to_root(out_dfs, config, central_cutset)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get BDT systematic')
    parser.add_argument('config', type=str, help='Path to the configuration file')
    parser.add_argument('--draw-only', action='store_true', help='Only draw the results')
    args = parser.parse_args()

    zfit.run.set_cpus_explicit(intra=20, inter=20)

    cut_variation(args.config, args.draw_only)

"""
Script for the raw-yield extraction of B mesons
"""

import argparse
import os

import numpy as np
import pandas as pd
import uproot
import yaml
from hist import Hist
from matplotlib.offsetbox import AnchoredText
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter

def create_hist(pt_lims, contents, errors, label_pt=r"$p_\mathrm{T}~(\mathrm{GeV}/c)$"):
    """
    Helper method to create histogram

    Parameters
    ----------

    - pt_lims (list): list of pt limits
    - contents (list): histogram contents
    - errors (list): histogram errors
    - label_pt (str): label for x axis

    Returns
    ----------
    - histogram (hist.Hist)

    """

    pt_cent = [(pt_min+pt_max)/2 for pt_min, pt_max in zip(pt_lims[:-1], pt_lims[1:])]
    histo = Hist.new.Var(pt_lims, name="x", label=label_pt).Weight()
    histo.fill(pt_cent, weight=contents)
    histo.view(flow=False).variance = np.array(errors)**2

    return histo


# pylint: disable=too-many-arguments
def add_info_on_canvas(axs, loc, system, pt_min, pt_max, fitter=None):
    """
    Helper method to add text on flarefly mass fit plot

    Parameters
    ----------
    - axs: matplotlib.figure.Axis
        Axis instance of the mass fit figure

    - loc: str
        Location of the info on the figure

    - system: str
        System (pp, MC pp)

    - pt_min: float
        Minimum pT value in the pT range

    - pt_max: float
        Maximum pT value in the pT range

    - fitter: F2MassFitter
        Fitter instance allowing to access chi2 and ndf if wanted
    """
    xspace = " "
    text = xspace
    if fitter is not None:
        chi2 = fitter.get_chi2()
        ndf = fitter.get_ndf()
        text += fr"$\chi^2 / \mathrm{{ndf}} =${chi2:.2f} / {ndf} $\simeq$ {chi2/ndf:.2f}""\n"

    text += "\n\n"
    text += xspace + system + ", " + r"$\sqrt{s} = 13.6$ TeV" + "\n"
    text += xspace + fr"{pt_min:.0f} < $p_{{\mathrm{{T}}}}$ < {pt_max:.0f} GeV/$c$, $|y|$ < 0.5""\n"

    anchored_text = AnchoredText(text, loc=loc, frameon=False)
    axs.add_artist(anchored_text)


def fit(config_file): # pylint: disable=too-many-locals,too-many-statements, too-many-branches
    """
    Main function for fitting

    Parameters
    ----------

    - config_file (string): config file name
    """

    with open(config_file, "r") as yml_cfg:  # pylint: disable=unspecified-encoding
        cfg = yaml.load(yml_cfg, yaml.FullLoader)

    with open(cfg["cutset_file_name"], "r") as yml_cfg:  # pylint: disable=unspecified-encoding
        cut_set = yaml.load(yml_cfg, yaml.FullLoader)

    particle = cfg["particle"]
    pdg_id = -1
    decay_channel = ""
    particle_name = ""
    flag_mc_match_rec = -999

    if particle == "Bplus":
        pdg_id = 521
        decay_channel = r"\overline{D}^{0} \pi^{\plus}"
        particle_name = "B^{+}"
        flag_mc_match_rec = 4 # prd = partly reco decays
    if particle == "B0":
        pdg_id = 511
        decay_channel = r"D^{-} \pi^{\plus}"
        particle_name = "B^{0}"
        flag_mc_match_rec = 8 # prd = partly reco decays

    pt_mins = cut_set["pt"]["mins"]
    pt_maxs = cut_set["pt"]["maxs"]
    pt_lims = pt_mins.copy()
    pt_lims.append(pt_maxs[-1])
    bdt_cut_mins = cut_set["ML_output"]["mins"]
    bdt_cut_maxs = cut_set["ML_output"]["maxs"]
    selection_string = ""
    for ipt, (pt_min, pt_max, bdt_cut_min, bdt_cut_max) in enumerate(zip(pt_mins, pt_maxs, bdt_cut_mins, bdt_cut_maxs)):
        if ipt == 0:
            selection_string += f"({pt_min} < fPt < {pt_max} and {bdt_cut_min} < ML_output < {bdt_cut_max})"
        else:
            selection_string += f" or ({pt_min} < fPt < {pt_max} and {bdt_cut_min} < ML_output < {bdt_cut_max})"

    # load data
    df = pd.DataFrame()
    for file in cfg["inputs"]["data"]:
        df = pd.concat([df, pd.read_parquet(file)])
    df.query(selection_string, inplace=True)

    # load mc
    df_mc = pd.DataFrame()
    for file in cfg["inputs"]["mc"]:
        df_mc = pd.concat([df_mc, pd.read_parquet(file)])
    df_mc.query(selection_string, inplace=True)
    df_mc_sig = df_mc.query("fFlagMcMatchRec == -1 or fFlagMcMatchRec == 1")

    dfs_prd_bkg, dfs_prd_bkg_orig, fracs_ptint = [], [], []
    correlated_bkgs = cfg["fit_configs"]["correlated_bkgs"]
    den_norm = len(df_mc_sig) * cfg["fit_configs"]["signal_br"]["pdg"] / cfg["fit_configs"]["signal_br"]["sim"]
    for bkg in correlated_bkgs:
        df_prd_bkg = df_mc.query(f"fFlagMcMatchRec == {flag_mc_match_rec} and "
                                 f"fPdgCodeBeautyMother == {bkg['beauty_id']} and "
                                 f"fPdgCodeCharmMother == {bkg['charm_id']}")
        dfs_prd_bkg_orig.append(df_prd_bkg)
        fracs_ptint.append(len(df_prd_bkg) * bkg["br_pdg"] / bkg["br_sim"] / den_norm)
    sum_fracs = sum(fracs_ptint)
    fracs_ptint_norm = [frac / sum_fracs for frac in fracs_ptint]

    dfs_prd_bkg_sampled = []
    if cfg["fit_configs"]["pt_int"]["bkg_templ_opt"] == 1:
        for frac, df_bkg in zip(fracs_ptint_norm, dfs_prd_bkg_orig):
            dfs_prd_bkg_sampled.append(df_bkg.sample(frac=frac, random_state=42))

        dfs_prd_bkg.append(pd.concat(dfs_prd_bkg_sampled))
    else:
        dfs_prd_bkg = dfs_prd_bkg_orig

    # define output file
    outdir = cfg["outputs"]["directory"]
    outfile_name = os.path.join(outdir,
                                f"{particle}_mass{cfg['outputs']['suffix']}.root")
    file_root = uproot.recreate(outfile_name)
    file_root.close()

    # we first perform a pT-integrated fit
    if cfg["fit_configs"]["pt_int"]["activate"]: #pylint:disable=too-many-nested-blocks
        # we first fit MC only
        data_hdl_mc = DataHandler(df_mc_sig, var_name="fM",
                                  limits=cfg["fit_configs"]["pt_int"]["mass_limits"],
                                  nbins=cfg["plot_style"]["pt_int"]["n_bins"])
        fitter_mc_ptint = F2MassFitter(data_hdl_mc,
                                       ["doublegaus"],
                                       ["nobkg"],
                                       name=f"{particle}_mc_ptint",
                                       label_signal_pdf=[rf"$\mathrm{{{particle_name}}}$ signal"])
        fitter_mc_ptint.set_signal_initpar(0, "sigma1", 0.03, limits=[0.01, 0.10])
        fitter_mc_ptint.set_signal_initpar(0, "sigma2", 0.085, limits=[0.01, 0.25])
        fitter_mc_ptint.set_particle_mass(0, pdg_id=pdg_id)
        result = fitter_mc_ptint.mass_zfit()
        if result.converged:
            fig, axs = fitter_mc_ptint.plot_mass_fit(style="ATLAS",
                                                     figsize=(8, 8),
                                                     axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)")
            add_info_on_canvas(axs, "upper left", "MC pp", pt_mins[0], pt_maxs[-1], fitter_mc_ptint)

            fig_res = fitter_mc_ptint.plot_raw_residuals(style="ATLAS",
                                                         figsize=(8, 8),
                                                         axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)")

            fig.savefig(os.path.join(outdir, f"{particle}_mass_ptint_MC.pdf"))
            fig_res.savefig(os.path.join(outdir, f"{particle}_massres_ptint_MC.pdf"))

        # then we fit data
        data_hdl = DataHandler(df, var_name="fM",
                               limits=cfg["fit_configs"]["pt_int"]["mass_limits"],
                               nbins=cfg["plot_style"]["pt_int"]["n_bins"])
        bkg_funcs = cfg["fit_configs"]["pt_int"]["bkg_funcs"]
        label_bkg_pdf = ["Comb. bkg"]
        data_hdls_prd_bkg = []
        if cfg["fit_configs"]["pt_int"]["use_bkg_templ"]:
            for i_bkg, (bkg, df_prd_bkg) in enumerate(zip(correlated_bkgs, dfs_prd_bkg)):
                data_hdls_prd_bkg.append(DataHandler(df_prd_bkg, var_name="fM",
                                                     limits=cfg["fit_configs"]["pt_int"]["mass_limits"],
                                                     nbins=cfg["plot_style"]["pt_int"]["n_bins"]))
                bkg_funcs.insert(i_bkg, "kde_grid")
                if cfg["fit_configs"]["pt_int"]["bkg_templ_opt"] == 0:
                    label_bkg_pdf.insert(i_bkg, bkg["name"])
                else:
                    label_bkg_pdf.insert(i_bkg, "Correlated backgrounds")

        fitter_ptint = F2MassFitter(data_hdl,
                                    cfg["fit_configs"]["pt_int"]["signal_funcs"],
                                    bkg_funcs,
                                    name=f"{particle}_ptint",
                                    label_signal_pdf=[rf"$\mathrm{{{particle_name}}}$ signal"],
                                    label_bkg_pdf=label_bkg_pdf)

        if cfg["fit_configs"]["pt_int"]["use_bkg_templ"]:
            for i_bkg, (bkg, data_hdl_prd_bkg) in enumerate(zip(correlated_bkgs, data_hdls_prd_bkg)):
                fitter_ptint.set_background_kde(i_bkg, data_hdl_prd_bkg)
                if i_bkg == 0:
                    if cfg["fit_configs"]["pt_int"]["fix_correlated_bkg_to_signal"]:
                        if cfg["fit_configs"]["pt_int"]["bkg_templ_opt"] == 0:
                            fitter_ptint.fix_bkg_frac_to_signal_pdf(i_bkg, 0, fracs_ptint[i_bkg])
                        else:
                            fitter_ptint.fix_bkg_frac_to_signal_pdf(i_bkg, 0, sum(fracs_ptint))
                    continue
                denom = (data_hdls_prd_bkg[0].get_norm() * correlated_bkgs[0]["br_pdg"] / correlated_bkgs[0]["br_sim"])
                fitter_ptint.fix_bkg_frac_to_bkg_pdf(
                    i_bkg, 0,
                    data_hdl_prd_bkg.get_norm() * bkg["br_pdg"] / bkg["br_sim"] / denom
                )

        fitter_ptint.set_signal_initpar(0, "sigma", 0.05, limits=[0.03, 0.06])
        fitter_ptint.set_particle_mass(0, pdg_id=pdg_id)
        icombbkg = len(dfs_prd_bkg)
        fitter_ptint.set_background_initpar(icombbkg, "c1", -0.05, limits=[-0.2, 0.])
        fitter_ptint.set_background_initpar(icombbkg, "c2", 0.008, limits=[0.000, 0.030])
        fitter_ptint.set_background_initpar(icombbkg, "lam", -1.2, limits=[-10., 10.])
        fitter_ptint.set_signal_initpar(0, "frac", 0.05, limits=[0., 1.])
        result = fitter_ptint.mass_zfit()
        if result.converged:
            fig, axs = fitter_ptint.plot_mass_fit(style="ATLAS",
                                                  figsize=(8, 8),
                                                  axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)",
                                                  show_extra_info=True,
                                                  extra_info_loc=["lower right", "lower left"])
            add_info_on_canvas(axs, "upper left", "pp", pt_mins[0], pt_maxs[-1])

            fig_res = fitter_ptint.plot_raw_residuals(style="ATLAS",
                                                      figsize=(8, 8),
                                                      axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)")

            fig.savefig(os.path.join(outdir, f"{particle}_mass_ptint.pdf"))
            fig_res.savefig(os.path.join(outdir, f"{particle}_massres_ptint.pdf"))

            fitter_ptint.dump_to_root(
                outfile_name, option="update", suffix="_ptint")

    raw_yields, raw_yields_unc = [], []
    signif, signif_unc, s_over_b, s_over_b_unc = [], [], [], []
    means, means_unc, sigmas, sigmas_unc = [], [], [], []
    means_mc, means_mc_unc, sigmas_mc, sigmas_mc_unc = [], [], [], []

    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)): #pylint:disable=too-many-nested-blocks

        # we first fit MC only
        df_mc_sig_pt = df_mc_sig.query(f"{pt_min} < fPt < {pt_max}")
        data_hdl_mc = DataHandler(df_mc_sig_pt, var_name="fM",
                                  limits=cfg["fit_configs"]["mass_limits"][ipt],
                                  nbins=cfg["plot_style"]["n_bins"][ipt])
        fitter_mc_pt = F2MassFitter(data_hdl_mc,
                                    ["doublegaus"],
                                    ["nobkg"],
                                    name=f"{particle}_mc_pt{pt_min:.0f}_{pt_max:.0f}",
                                    label_signal_pdf=[rf"$\mathrm{{{particle_name}}}$ signal"])
        fitter_mc_pt.set_signal_initpar(0, "sigma1", 0.03, limits=[0.01, 0.10])
        fitter_mc_pt.set_signal_initpar(0, "sigma2", 0.085, limits=[0.01, 0.25])
        fitter_mc_pt.set_particle_mass(0, pdg_id=pdg_id)
        result = fitter_mc_pt.mass_zfit()
        if result.converged:
            fig, axs = fitter_mc_pt.plot_mass_fit(style="ATLAS",
                                                  figsize=(8, 8),
                                                  axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)")
            add_info_on_canvas(axs, "upper left", "MC pp", pt_min, pt_max, fitter_mc_pt)

            fig_res = fitter_mc_pt.plot_raw_residuals(style="ATLAS",
                                                      figsize=(8, 8),
                                                      axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)")

            fig.savefig(os.path.join(outdir, f"{particle}_mass_pt{pt_min:.0f}_{pt_max:.0f}_MC.pdf"))
            fig_res.savefig(os.path.join(outdir, f"{particle}_massres_pt{pt_min:.0f}_{pt_max:.0f}_MC.pdf"))

            mean, mean_unc = fitter_mc_pt.get_signal_parameter(0, "mu")
            sigma1, sigma1_unc = fitter_mc_pt.get_signal_parameter(0, "sigma1")
            sigma2, sigma2_unc = fitter_mc_pt.get_signal_parameter(0, "sigma2")
            if sigma1 < sigma2:
                sigma = sigma1
                sigma_unc = sigma1_unc
            else:
                sigma = sigma2
                sigma_unc = sigma2_unc

            means_mc.append(mean)
            means_mc_unc.append(mean_unc)
            sigmas_mc.append(sigma)
            sigmas_mc_unc.append(sigma_unc)

        # then we fit data
        df_pt = df.query(f"{pt_min} < fPt < {pt_max}")
        data_hdl = DataHandler(df_pt, var_name="fM",
                               limits=cfg["fit_configs"]["mass_limits"][ipt],
                               nbins=cfg["plot_style"]["n_bins"][ipt])

        den_norm = data_hdl_mc.get_norm() * cfg["fit_configs"]["signal_br"]["pdg"] / \
            cfg["fit_configs"]["signal_br"]["sim"]

        bkg_funcs = cfg["fit_configs"]["bkg_funcs"][ipt]
        label_bkg_pdf = ["Comb. bkg"]
        data_hdls_prd_bkg = []
        dfs_prd_bkg_orig_pt = [df.query(f"{pt_min} < fPt < {pt_max}") for df in dfs_prd_bkg_orig]
        fracs_pt = []
        for bkg, df_prd_bkg_orig_pt in zip(correlated_bkgs, dfs_prd_bkg_orig_pt):
            fracs_pt.append(len(df_prd_bkg_orig_pt) * bkg["br_pdg"] / bkg["br_sim"] / den_norm)
        sum_fracs = sum(fracs_pt)
        fracs_pt_norm = [frac / sum_fracs for frac in fracs_pt]

        dfs_prd_bkg_pt = []
        dfs_prd_bkg_sampled_pt = []
        if cfg["fit_configs"]["pt_int"]["bkg_templ_opt"] == 1:
            for frac, df_bkg in zip(fracs_pt_norm, dfs_prd_bkg_orig_pt):
                dfs_prd_bkg_sampled_pt.append(df_bkg.sample(frac=frac, random_state=42))

            dfs_prd_bkg_pt.append(pd.concat(dfs_prd_bkg_sampled_pt))
        else:
            dfs_prd_bkg_pt = dfs_prd_bkg_orig_pt

        for i_bkg, (bkg, df_prd_bkg) in enumerate(zip(correlated_bkgs, dfs_prd_bkg_pt)):
            data_hdls_prd_bkg.append(DataHandler(df_prd_bkg, var_name="fM",
                                                    limits=cfg["fit_configs"]["pt_int"]["mass_limits"],
                                                    nbins=cfg["plot_style"]["pt_int"]["n_bins"]))
            bkg_funcs.insert(i_bkg, "kde_grid")
            if cfg["fit_configs"]["bkg_templ_opt"][ipt] == 0:
                label_bkg_pdf.insert(i_bkg, bkg["name"])
            else:
                label_bkg_pdf.insert(i_bkg, "Correlated backgrounds")                   

        fitter_pt = F2MassFitter(data_hdl,
                                 cfg["fit_configs"]["signal_funcs"][ipt],
                                 bkg_funcs,
                                 name=f"{particle}_pt{pt_min:.0f}_{pt_max:.0f}",
                                 label_signal_pdf=[rf"$\mathrm{{{particle_name}}}$ signal"],
                                 label_bkg_pdf=label_bkg_pdf
                                 )
        if cfg["fit_configs"]["use_bkg_templ"][ipt]:
            for i_bkg, (bkg, data_hdl_prd_bkg) in enumerate(zip(correlated_bkgs, data_hdls_prd_bkg)):
                fitter_pt.set_background_kde(i_bkg, data_hdl_prd_bkg)
                if i_bkg == 0:
                    if cfg["fit_configs"]["fix_correlated_bkg_to_signal"][ipt]:
                        if cfg["fit_configs"]["bkg_templ_opt"][ipt] == 0:
                            fitter_pt.fix_bkg_frac_to_signal_pdf(i_bkg, 0, fracs_pt[i_bkg])
                        else:
                            fitter_pt.fix_bkg_frac_to_signal_pdf(i_bkg, 0, sum(fracs_pt))
                    continue
                denom = (data_hdls_prd_bkg[0].get_norm() * correlated_bkgs[0]["br_pdg"] / correlated_bkgs[0]["br_sim"])
                fitter_pt.fix_bkg_frac_to_bkg_pdf(
                    i_bkg, 0,
                    data_hdl_prd_bkg.get_norm() * bkg["br_pdg"] / bkg["br_sim"] / denom
                )

        fitter_pt.set_signal_initpar(0, "sigma", sigmas_mc[ipt], limits=[0.01, 0.1])
        fitter_pt.set_particle_mass(0, pdg_id=pdg_id, limits=[5., 5.56])
        fitter_pt.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])
        if not cfg["fit_configs"]["fix_correlated_bkg_to_signal"][ipt]:
            fitter_pt.set_background_initpar(0, "frac", 0.05, limits=[0., 1.])
        icombbkg = len(dfs_prd_bkg_pt)
        fitter_pt.set_background_initpar(icombbkg, "lam", -1.2, limits=[-10., 10.])
        fitter_pt.set_background_initpar(icombbkg, "c1", -0.05, limits=[-0.2, 0.])
        fitter_pt.set_background_initpar(icombbkg, "c2", 0.008, limits=[0.000, 0.03])
        result = fitter_pt.mass_zfit()
        if result.converged:
            fig, axs = fitter_pt.plot_mass_fit(style="ATLAS",
                                               figsize=(8, 8),
                                               axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)",
                                               show_extra_info=True,
                                               extra_info_loc=["lower right", "lower left"])
            add_info_on_canvas(axs, "upper left", "pp", pt_min, pt_max)

            fig_res = fitter_pt.plot_raw_residuals(style="ATLAS",
                                                   figsize=(8, 8),
                                                   axis_title=rf"$M(\mathrm{{{decay_channel}}})$ (GeV/$c^2$)")

            fig.savefig(os.path.join(outdir, f"{particle}_mass_pt{pt_min:.0f}_{pt_max:.0f}.pdf"))
            fig_res.savefig(os.path.join(outdir, f"{particle}_massres_pt{pt_min:.0f}_{pt_max:.0f}.pdf"))

            rawy, rawy_unc = fitter_pt.get_raw_yield(0)
            sign, sign_unc = fitter_pt.get_significance(0)
            soverb, soverb_unc = fitter_pt.get_signal_over_background(0)
            mean, mean_unc = fitter_pt.get_signal_parameter(0, "mu")
            sigma, sigma_unc = fitter_pt.get_signal_parameter(0, "sigma")

            raw_yields.append(rawy)
            raw_yields_unc.append(rawy_unc)
            signif.append(sign)
            signif_unc.append(sign_unc)
            s_over_b.append(soverb)
            s_over_b_unc.append(soverb_unc)
            means.append(mean)
            means_unc.append(mean_unc)
            sigmas.append(sigma)
            sigmas_unc.append(sigma_unc)

            fitter_pt.dump_to_root(
                outfile_name, option="update", suffix=f"_pt{pt_min:.0f}_{pt_max:.0f}")

    file_root = uproot.update(outfile_name)
    file_root["h_rawyields"] = create_hist(pt_lims, raw_yields, raw_yields_unc)
    file_root["h_significance"] = create_hist(pt_lims, signif, signif_unc)
    file_root["h_soverb"] = create_hist(pt_lims, s_over_b, s_over_b_unc)
    file_root["h_means"] = create_hist(pt_lims, means, means_unc)
    file_root["h_sigmas"] = create_hist(pt_lims, sigmas, sigmas_unc)
    file_root["h_means_mc"] = create_hist(pt_lims, means_mc, means_mc_unc)
    file_root["h_sigmas_mc"] = create_hist(pt_lims, sigmas_mc, sigmas_mc_unc)
    file_root.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-c", metavar="text", default="config_fit.yml",
                        help="yaml config file for fit", required=True)
    args = parser.parse_args()
    fit(args.config)

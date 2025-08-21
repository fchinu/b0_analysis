"""
Script to plot FONLL predictions for B mesons
"""

import sys
import argparse
import math
import yaml
import pandas as pd
import ROOT

sys.path.append('../utils') # pylint: disable=wrong-import-position
from style_formatter import set_global_style, set_object_style # pylint: disable=import-error


# pylint: disable=too-many-locals
def rebin_df(df, rebin):
    """
    Function to rebin a pandas dataframe
    """

    ptmins = rebin.copy()
    ptmaxs = rebin.copy()
    ptmins.pop()
    ptmaxs.pop(0)

    #first check pt consistency
    ptmins_orig = df["ptmin"].to_numpy()
    ptmaxs_orig = df["ptmax"].to_numpy()
    ptmins_ok = all(any(
        math.isclose(pt, pt_orig, rel_tol=1e-6) for pt_orig in ptmins_orig) for pt in ptmins)
    ptmaxs_ok = all(any(
        math.isclose(pt, pt_orig, rel_tol=1e-6) for pt_orig in ptmaxs_orig) for pt in ptmaxs)

    if not ptmins_ok or not ptmaxs_ok:
        print("ERROR: pt rebin from config not consistent with pt bins in input files")
        sys.exit()

    n_ptbins = len(ptmins)

    l_xsec, l_xsec_min, l_xsec_max,  = ([0. for _ in range(n_ptbins)] for _ in range(3))
    l_min_sc, l_max_sc, l_min_mass = ([0. for _ in range(n_ptbins)] for _ in range(3))
    l_max_mass, l_min_pdf, l_max_pdf = ([0. for _ in range(n_ptbins)] for _ in range(3))
    l_fr_dot5_dot5, l_fr_2_2, l_fr_2_1 = ([0. for _ in range(n_ptbins)] for _ in range(3))
    l_fr_1_2, l_fr_1_dot5, l_fr_dot5_1 = ([0. for _ in range(n_ptbins)] for _ in range(3))

    for ipt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
        for (ptmin_orig, ptmax_orig, xsec, xsec_min, xsec_max, min_sc, max_sc,
             min_mass, max_mass, min_pdf, max_pdf, fr_dot5_dot5, fr_2_2, fr_2_1, fr_1_2,
             fr_1_dot5, fr_dot5_1) in zip(df["ptmin"].to_numpy(),
                                          df["ptmax"].to_numpy(),
                                          df["central"].to_numpy(),
                                          df["min"].to_numpy(),
                                          df["max"].to_numpy(),
                                          df["min_sc"].to_numpy(),
                                          df["max_sc"].to_numpy(),
                                          df["min_mass"].to_numpy(),
                                          df["max_mass"].to_numpy(),
                                          df["min_pdf"].to_numpy(),
                                          df["max_pdf"].to_numpy(),
                                          df["fr_dot5_dot5"].to_numpy(),
                                          df["fr_2_2"].to_numpy(),
                                          df["fr_2_1"].to_numpy(),
                                          df["fr_1_2"].to_numpy(),
                                          df["fr_1_dot5"].to_numpy(),
                                          df["fr_dot5_1"].to_numpy()):

            if (ptmin_orig > ptmin or math.isclose(ptmin_orig, ptmin, rel_tol=1e-6)) and \
                (ptmax_orig < ptmax or math.isclose(ptmax_orig, ptmax, rel_tol=1e-6)):
                l_xsec[ipt] += xsec
                l_xsec_min[ipt] += xsec_min
                l_xsec_max[ipt] += xsec_max
                l_min_sc[ipt] += min_sc
                l_max_sc[ipt] += max_sc
                l_min_mass[ipt] += min_mass
                l_max_mass[ipt] += max_mass
                l_min_pdf[ipt] += min_pdf
                l_max_pdf[ipt] += max_pdf
                l_fr_dot5_dot5[ipt] += fr_dot5_dot5
                l_fr_2_2[ipt] += fr_2_2
                l_fr_2_1[ipt] += fr_2_1
                l_fr_1_2[ipt] += fr_1_2
                l_fr_1_dot5[ipt] += fr_1_dot5
                l_fr_dot5_1[ipt] += fr_dot5_1

    df_rebin = pd.DataFrame({
        "ptmin": ptmins,
        "ptmax": ptmaxs,
        "central": l_xsec,
        "min": l_xsec_min,
        "max": l_xsec_max,
        "min_sc": l_min_sc,
        "max_sc": l_max_sc,
        "min_mass": l_min_mass,
        "max_mass": l_max_mass,
        "min_pdf": l_min_pdf,
        "max_pdf": l_max_pdf,
        "fr_dot5_dot5": l_fr_dot5_dot5,
        "fr_2_2": l_fr_2_2,
        "fr_2_1": l_fr_2_1,
        "fr_1_2": l_fr_1_2,
        "fr_1_dot5": l_fr_1_dot5,
        "fr_dot5_1": l_fr_dot5_1
    })

    return df_rebin


def get_rapidity_interval_and_ff(file_name):
    """
    Function to retrieve rapidity interval and fragmentation fraction set in the FONLL website
    """

    y_min, y_max, ff = -999., -999., 0.
    with open(file_name, "r") as file:
        for line in file:
            line_stripped = line.strip()
            if "ymin" in line_stripped:
                y_min = float(line_stripped.split(sep=" ")[-1])
            elif "ymax" in line_stripped:
                y_max = float(line_stripped.split(sep=" ")[-1])
            elif "BR(q->meson)" in line_stripped:
                ff = float(line_stripped.split(sep=" ")[-1])

    return [y_min, y_max], ff


def convert_to_graph(df, deltay, ff, graph_name, graph_color):
    """
    Simple function to convert a pandas dataframe to a graph
    """

    graph = ROOT.TGraphAsymmErrors(len(df))
    graph.SetNameTitle(graph_name, ";#it{p}_{T} (GeV/#it{c});"
                       "d^{2}#sigma/d#it{p}_{T}d#it{y} (pb #it{c}/GeV)")
    set_object_style(graph, color=graph_color, alpha=0.5)

    for ipt, (ptmin, ptmax, xsec, xsecmin, xsecmax) in enumerate(zip(
        df["ptmin"].to_numpy(),
        df["ptmax"].to_numpy(),
        df["central"].to_numpy(),
        df["min"].to_numpy(),
        df["max"].to_numpy())):

        deltapt = ptmax - ptmin
        ptcent = (ptmax + ptmin) / 2
        ptunc = deltapt / 2
        diffxsec = xsec / deltapt / deltay * ff
        diffxsec_min = xsecmin / deltapt / deltay * ff
        diffxsec_max = xsecmax / deltapt / deltay * ff

        graph.SetPoint(ipt, ptcent, diffxsec)
        graph.SetPointError(ipt, ptunc, ptunc, diffxsec-diffxsec_min, diffxsec_max-diffxsec)

    return graph


# pylint: disable=too-many-arguments
def get_ratio_fwd_mid(df_mid, df_fwd, deltay_mid, deltay_fwd, ff_mid, ff_fwd, graph_name, graph_color):
    """
    Function to compute mid / fwd FONLL ratio
    """

    graph_ratio = ROOT.TGraphAsymmErrors(len(df_mid))
    graph_ratio.SetNameTitle(graph_name,
                             ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (mid/fwd)")
    set_object_style(graph_ratio, color=graph_color, alpha=0.5)
    for ipt, (ptmin, ptmax, cent_mid, min_sc_mid, max_sc_mid,
              min_mass_mid, max_mass_mid, min_pdf_mid, max_pdf_mid,
              fr_dot5_dot5_mid, fr_2_2_mid, fr_2_1_mid, fr_1_2_mid,
              fr_1_dot5_mid, fr_dot5_1_mid, cent_fwd, min_sc_fwd, max_sc_fwd,
              min_mass_fwd, max_mass_fwd, min_pdf_fwd, max_pdf_fwd,
              fr_dot5_dot5_fwd, fr_2_2_fwd, fr_2_1_fwd, fr_1_2_fwd,
              fr_1_dot5_fwd, fr_dot5_1_fwd) in enumerate(zip(df_mid["ptmin"].to_numpy(),
                                                             df_mid["ptmax"].to_numpy(),
                                                             df_mid["central"].to_numpy(),
                                                             df_mid["min_sc"].to_numpy(),
                                                             df_mid["max_sc"].to_numpy(),
                                                             df_mid["min_mass"].to_numpy(),
                                                             df_mid["max_mass"].to_numpy(),
                                                             df_mid["min_pdf"].to_numpy(),
                                                             df_mid["max_pdf"].to_numpy(),
                                                             df_mid["fr_dot5_dot5"].to_numpy(),
                                                             df_mid["fr_2_2"].to_numpy(),
                                                             df_mid["fr_2_1"].to_numpy(),
                                                             df_mid["fr_1_2"].to_numpy(),
                                                             df_mid["fr_1_dot5"].to_numpy(),
                                                             df_mid["fr_dot5_1"].to_numpy(),
                                                             df_fwd["central"].to_numpy(),
                                                             df_fwd["min_sc"].to_numpy(),
                                                             df_fwd["max_sc"].to_numpy(),
                                                             df_fwd["min_mass"].to_numpy(),
                                                             df_fwd["max_mass"].to_numpy(),
                                                             df_fwd["min_pdf"].to_numpy(),
                                                             df_fwd["max_pdf"].to_numpy(),
                                                             df_fwd["fr_dot5_dot5"].to_numpy(),
                                                             df_fwd["fr_2_2"].to_numpy(),
                                                             df_fwd["fr_2_1"].to_numpy(),
                                                             df_fwd["fr_1_2"].to_numpy(),
                                                             df_fwd["fr_1_dot5"].to_numpy(),
                                                             df_fwd["fr_dot5_1"].to_numpy())):
        ptcent = (ptmax + ptmin) / 2
        ptunc = (ptmax - ptmin) / 2
        norm = 1. / deltay_mid * deltay_fwd * ff_mid / ff_fwd
        ratio = cent_mid / cent_fwd * norm
        ratios_vars = [
            min_sc_mid / min_sc_fwd * norm,
            max_sc_mid / max_sc_fwd * norm,
            min_mass_mid / min_mass_fwd * norm,
            max_mass_mid / max_mass_fwd * norm,
            min_pdf_mid / min_pdf_fwd * norm,
            max_pdf_mid / max_pdf_fwd * norm,
            fr_dot5_dot5_mid / fr_dot5_dot5_fwd * norm,
            fr_2_2_mid / fr_2_2_fwd * norm,
            fr_2_1_mid / fr_2_1_fwd * norm,
            fr_1_2_mid / fr_1_2_fwd * norm,
            fr_1_dot5_mid / fr_1_dot5_fwd * norm,
            fr_dot5_1_mid / fr_dot5_1_fwd * norm,
        ]

        graph_ratio.SetPoint(ipt, ptcent, ratio)
        graph_ratio.SetPointError(ipt, ptunc, ptunc,
                                  ratio - min(ratios_vars), max(ratios_vars) - ratio)

    return graph_ratio


# pylint: disable=too-many-arguments
def get_double_ratio_fwd_mid(df_mid_c, df_fwd_c, df_mid_b, df_fwd_b, deltay_mid_c, deltay_fwd_c, deltay_mid_b, deltay_fwd_b, ff_mid_c, ff_fwd_c, ff_mid_b, ff_fwd_b, graph_name, graph_color):
    """
    Function to compute mid / fwd FONLL ratio
    """

    graph_ratio = ROOT.TGraphAsymmErrors(len(df_mid_c))
    graph_ratio.SetNameTitle(graph_name,
                             ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (mid/fwd)")
    set_object_style(graph_ratio, color=graph_color, alpha=0.5)
    for ipt, (ptmin, ptmax, cent_mid_c, min_sc_mid_c, max_sc_mid_c,
              min_mass_mid_c, max_mass_mid_c, min_pdf_mid_c, max_pdf_mid_c,
              fr_dot5_dot5_mid_c, fr_2_2_mid_c, fr_2_1_mid_c, fr_1_2_mid_c,
              fr_1_dot5_mid_c, fr_dot5_1_mid_c, cent_fwd_c, min_sc_fwd_c, max_sc_fwd_c,
              min_mass_fwd_c, max_mass_fwd_c, min_pdf_fwd_c, max_pdf_fwd_c,
              fr_dot5_dot5_fwd_c, fr_2_2_fwd_c, fr_2_1_fwd_c, fr_1_2_fwd_c,
              fr_1_dot5_fwd_c, fr_dot5_1_fwd_c, cent_mid_b, min_sc_mid_b, max_sc_mid_b,
              min_mass_mid_b, max_mass_mid_b, min_pdf_mid_b, max_pdf_mid_b,
              fr_dot5_dot5_mid_b, fr_2_2_mid_b, fr_2_1_mid_b, fr_1_2_mid_b,
              fr_1_dot5_mid_b, fr_dot5_1_mid_b, cent_fwd_b, min_sc_fwd_b, max_sc_fwd_b,
              min_mass_fwd_b, max_mass_fwd_b, min_pdf_fwd_b, max_pdf_fwd_b,
              fr_dot5_dot5_fwd_b, fr_2_2_fwd_b, fr_2_1_fwd_b, fr_1_2_fwd_b,
              fr_1_dot5_fwd_b, fr_dot5_1_fwd_b) in enumerate(zip(df_mid_c["ptmin"].to_numpy(),
                                                                 df_mid_c["ptmax"].to_numpy(),
                                                                 df_mid_c["central"].to_numpy(),
                                                                 df_mid_c["min_sc"].to_numpy(),
                                                                 df_mid_c["max_sc"].to_numpy(),
                                                                 df_mid_c["min_mass"].to_numpy(),
                                                                 df_mid_c["max_mass"].to_numpy(),
                                                                 df_mid_c["min_pdf"].to_numpy(),
                                                                 df_mid_c["max_pdf"].to_numpy(),
                                                                 df_mid_c["fr_dot5_dot5"].to_numpy(),
                                                                 df_mid_c["fr_2_2"].to_numpy(),
                                                                 df_mid_c["fr_2_1"].to_numpy(),
                                                                 df_mid_c["fr_1_2"].to_numpy(),
                                                                 df_mid_c["fr_1_dot5"].to_numpy(),
                                                                 df_mid_c["fr_dot5_1"].to_numpy(),
                                                                 df_fwd_c["central"].to_numpy(),
                                                                 df_fwd_c["min_sc"].to_numpy(),
                                                                 df_fwd_c["max_sc"].to_numpy(),
                                                                 df_fwd_c["min_mass"].to_numpy(),
                                                                 df_fwd_c["max_mass"].to_numpy(),
                                                                 df_fwd_c["min_pdf"].to_numpy(),
                                                                 df_fwd_c["max_pdf"].to_numpy(),
                                                                 df_fwd_c["fr_dot5_dot5"].to_numpy(),
                                                                 df_fwd_c["fr_2_2"].to_numpy(),
                                                                 df_fwd_c["fr_2_1"].to_numpy(),
                                                                 df_fwd_c["fr_1_2"].to_numpy(),
                                                                 df_fwd_c["fr_1_dot5"].to_numpy(),
                                                                 df_fwd_c["fr_dot5_1"].to_numpy(),
                                                                 df_mid_b["central"].to_numpy(),
                                                                 df_mid_b["min_sc"].to_numpy(),
                                                                 df_mid_b["max_sc"].to_numpy(),
                                                                 df_mid_b["min_mass"].to_numpy(),
                                                                 df_mid_b["max_mass"].to_numpy(),
                                                                 df_mid_b["min_pdf"].to_numpy(),
                                                                 df_mid_b["max_pdf"].to_numpy(),
                                                                 df_mid_b["fr_dot5_dot5"].to_numpy(),
                                                                 df_mid_b["fr_2_2"].to_numpy(),
                                                                 df_mid_b["fr_2_1"].to_numpy(),
                                                                 df_mid_b["fr_1_2"].to_numpy(),
                                                                 df_mid_b["fr_1_dot5"].to_numpy(),
                                                                 df_mid_b["fr_dot5_1"].to_numpy(),
                                                                 df_fwd_b["central"].to_numpy(),
                                                                 df_fwd_b["min_sc"].to_numpy(),
                                                                 df_fwd_b["max_sc"].to_numpy(),
                                                                 df_fwd_b["min_mass"].to_numpy(),
                                                                 df_fwd_b["max_mass"].to_numpy(),
                                                                 df_fwd_b["min_pdf"].to_numpy(),
                                                                 df_fwd_b["max_pdf"].to_numpy(),
                                                                 df_fwd_b["fr_dot5_dot5"].to_numpy(),
                                                                 df_fwd_b["fr_2_2"].to_numpy(),
                                                                 df_fwd_b["fr_2_1"].to_numpy(),
                                                                 df_fwd_b["fr_1_2"].to_numpy(),
                                                                 df_fwd_b["fr_1_dot5"].to_numpy(),
                                                                 df_fwd_b["fr_dot5_1"].to_numpy())):
        ptcent = (ptmax + ptmin) / 2
        ptunc = (ptmax - ptmin) / 2
        norm_b = 1. / deltay_mid_b * deltay_fwd_b * ff_mid_b / ff_fwd_b
        norm_c = 1. / deltay_mid_c * deltay_fwd_c * ff_mid_c / ff_fwd_c
        ratio = cent_mid_b / cent_fwd_b * norm_b / (cent_mid_c / cent_fwd_c * norm_c)
        ratios_vars = [
            min_sc_mid_b / min_sc_fwd_b * norm_b / (min_sc_mid_c / min_sc_fwd_c * norm_c),
            max_sc_mid_b / max_sc_fwd_b * norm_b / (max_sc_mid_c / max_sc_fwd_c * norm_c),
            min_mass_mid_b / min_mass_fwd_b * norm_b / (min_mass_mid_c / min_mass_fwd_c * norm_c),
            max_mass_mid_b / max_mass_fwd_b * norm_b / (max_mass_mid_c / max_mass_fwd_c * norm_c),
            min_pdf_mid_b / min_pdf_fwd_b * norm_b / (min_pdf_mid_c / min_pdf_fwd_c * norm_c),
            max_pdf_mid_b / max_pdf_fwd_b * norm_b / (max_pdf_mid_c / max_pdf_fwd_c * norm_c),
            fr_dot5_dot5_mid_b / fr_dot5_dot5_fwd_b * norm_b / (fr_dot5_dot5_mid_c / fr_dot5_dot5_fwd_c * norm_c),
            fr_2_2_mid_b / fr_2_2_fwd_b * norm_b / (fr_2_2_mid_c / fr_2_2_fwd_c * norm_c),
            fr_2_1_mid_b / fr_2_1_fwd_b * norm_b / (fr_2_1_mid_c / fr_2_1_fwd_c * norm_c),
            fr_1_2_mid_b / fr_1_2_fwd_b * norm_b / (fr_1_2_mid_c / fr_1_2_fwd_c * norm_c),
            fr_1_dot5_mid_b / fr_1_dot5_fwd_b * norm_b / (fr_1_dot5_mid_c / fr_1_dot5_fwd_c * norm_c),
            fr_dot5_1_mid_b / fr_dot5_1_fwd_b * norm_b / (fr_dot5_1_mid_c / fr_dot5_1_fwd_c * norm_c),
        ]

        graph_ratio.SetPoint(ipt, ptcent, ratio)
        graph_ratio.SetPointError(ipt, ptunc, ptunc,
                                  ratio - min(ratios_vars), max(ratios_vars) - ratio)

    return graph_ratio


def plot(config_file):
    """
    Main function to produce graphs and plots 
    """

    set_global_style(padleftmargin=0.18, padbottommargin=0.14)

    with open(config_file, "r") as yml_cfg:  # pylint: disable=unspecified-encoding
        cfg = yaml.load(yml_cfg, yaml.FullLoader)

    col_names = ["ptmin", "ptmax", "central", "min", "max", "min_sc", "max_sc",
                 "min_mass", "max_mass", "min_pdf", "max_pdf", "fr_dot5_dot5", "fr_2_2",
                 "fr_2_1", "fr_1_2", "fr_1_dot5", "fr_dot5_1"]
    df_mid_b = pd.read_csv(cfg["inputs"]["beauty"]["mid"], names=col_names, comment="#", sep=" ")
    df_fwd_b = pd.read_csv(cfg["inputs"]["beauty"]["fwd"], names=col_names, comment="#", sep=" ")
    df_mid_c = pd.read_csv(cfg["inputs"]["charm"]["mid"], names=col_names, comment="#", sep=" ")
    df_fwd_c = pd.read_csv(cfg["inputs"]["charm"]["fwd"], names=col_names, comment="#", sep=" ")

    if cfg["rebin_pt"]["beauty"]["enable"]:
        df_mid_b = rebin_df(df_mid_b, cfg["rebin_pt"]["beauty"]["ptlimits"])
        df_fwd_b = rebin_df(df_fwd_b, cfg["rebin_pt"]["beauty"]["ptlimits"])
    if cfg["rebin_pt"]["charm"]["enable"]:
        df_mid_c = rebin_df(df_mid_c, cfg["rebin_pt"]["charm"]["ptlimits"])
        df_fwd_c = rebin_df(df_fwd_c, cfg["rebin_pt"]["charm"]["ptlimits"])

    y_mid_b, ff_mid_b = get_rapidity_interval_and_ff(cfg["inputs"]["beauty"]["mid"])
    y_fwd_b, ff_fwd_b = get_rapidity_interval_and_ff(cfg["inputs"]["beauty"]["fwd"])
    y_mid_c, ff_mid_c = get_rapidity_interval_and_ff(cfg["inputs"]["charm"]["mid"])
    y_fwd_c, ff_fwd_c = get_rapidity_interval_and_ff(cfg["inputs"]["charm"]["fwd"])
    ff_mid_b = cfg["frag_fracs"]["beauty"]["mid"] / ff_mid_b
    ff_fwd_b = cfg["frag_fracs"]["beauty"]["fwd"] / ff_fwd_b
    ff_mid_c = cfg["frag_fracs"]["charm"]["mid"] / ff_mid_c
    ff_fwd_c = cfg["frag_fracs"]["charm"]["fwd"] / ff_fwd_c
    deltay_mid_b = y_mid_b[1] - y_mid_b[0]
    deltay_fwd_b = y_fwd_b[1] - y_fwd_b[0]
    deltay_mid_c = y_mid_c[1] - y_mid_c[0]
    deltay_fwd_c = y_fwd_c[1] - y_fwd_c[0]

    graph_mid_b = convert_to_graph(df_mid_b, deltay_mid_b, ff_mid_b, "graph_ratio_mid_b", ROOT.kAzure+4)
    graph_fwd_b = convert_to_graph(df_fwd_b, deltay_fwd_b, ff_fwd_b, "graph_ratio_fwd_b", ROOT.kBlue+3)
    graph_ratio_b = get_ratio_fwd_mid(df_mid_b, df_fwd_b, deltay_mid_b, deltay_fwd_b,
                                      ff_mid_b, ff_fwd_b, "graph_ratio_mid_fwd_b", ROOT.kAzure+4)
    graph_mid_c = convert_to_graph(df_mid_c, deltay_mid_c, ff_mid_c, "graph_ratio_mid_c", ROOT.kRed+2)
    graph_fwd_c = convert_to_graph(df_fwd_c, deltay_fwd_c, ff_fwd_c, "graph_ratio_fwd_c", ROOT.kRed+3)
    graph_ratio_c = get_ratio_fwd_mid(df_mid_c, df_fwd_c, deltay_mid_c, deltay_fwd_c,
                                      ff_mid_c, ff_fwd_c, "graph_ratio_mid_fwd_c", ROOT.kRed+2)

    graph_double_ratio = None
    if cfg["rebin_pt"]["beauty"]["ptlimits"] == cfg["rebin_pt"]["charm"]["ptlimits"]:
        graph_double_ratio = get_double_ratio_fwd_mid(df_mid_c, df_fwd_c, df_mid_b, df_fwd_b,
                                                      deltay_mid_c, deltay_fwd_c, deltay_mid_b, deltay_fwd_b,
                                                      ff_mid_c, ff_fwd_c, ff_mid_b, ff_fwd_b,
                                                      "graph_double_ratio_mid_fwd", ROOT.kAzure+4)

    line_at_one = ROOT.TLine(0., 1., graph_mid_b.GetPointX(graph_mid_b.GetN()-1) + 10, 1.)
    line_at_one.SetLineWidth(1)
    line_at_one.SetLineStyle(9)
    line_at_one.SetLineColor(ROOT.kGray + 1)

    leg = ROOT.TLegend(0.5, 0.7, 0.9, 0.85)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(graph_mid_b, cfg["legend"]["mid"], "fp")
    leg.AddEntry(graph_fwd_b, cfg["legend"]["fwd"], "fp")

    leg_ratio = ROOT.TLegend(0.3, 0.2, 0.9, 0.35)
    leg_ratio.SetTextSize(0.035)
    leg_ratio.SetFillStyle(0)
    leg_ratio.SetBorderSize(0)
    leg_ratio.AddEntry(graph_ratio_c, f"D {y_mid_c[0]:.1f} < #it{{y}} < {y_mid_c[1]:.1f} / {y_fwd_c[0]:.1f} < #it{{y}} < {y_fwd_c[1]:.1f}", "fp")
    leg_ratio.AddEntry(graph_ratio_b, f"B {y_mid_b[0]:.1f} < #it{{y}} < {y_mid_b[1]:.1f} / {y_fwd_b[0]:.1f} < #it{{y}} < {y_fwd_b[1]:.1f}", "fp")

    canvas = ROOT.TCanvas("canv_fonll", "", 600, 600)
    canvas.cd().DrawFrame(0., 1.e4,
                           graph_mid_b.GetPointX(graph_mid_b.GetN()-1)+10, 1.e8,
                           ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (pb #it{c}/GeV)")
    canvas.cd().SetLogy()
    graph_mid_b.Draw("2pz")
    graph_fwd_b.Draw("2pz")
    leg.Draw()
    canvas.SaveAs(cfg["output"].replace(".root", ".pdf"))

    canvas_ratio = ROOT.TCanvas("canvas_ratio", "", 600, 600)
    frame_ratio = canvas_ratio.cd().DrawFrame(
        0., 0.,
        graph_mid_b.GetPointX(graph_mid_b.GetN()-1)+10, 2.0,
        ";#it{p}_{T} (GeV/#it{c});Ratio mid / fwd"
    )
    line_at_one.Draw()
    frame_ratio.GetYaxis().SetDecimals()
    graph_ratio_b.Draw("2pz")
    graph_ratio_c.Draw("2pz")
    leg_ratio.Draw()
    canvas_ratio.SaveAs(cfg["output"].replace(".root", "_ratio.pdf"))

    if graph_double_ratio is not None:
        canvas_double_ratio = ROOT.TCanvas("canvas_double_ratio", "", 600, 600)
        frame_ratio = canvas_double_ratio.cd().DrawFrame(
            0., 0.,
            graph_mid_b.GetPointX(graph_mid_b.GetN()-1)+10, 2.0,
            ";#it{p}_{T} (GeV/#it{c});Ratio mid / fwd (beauty / charm)"
        )
        line_at_one.Draw()
        frame_ratio.GetYaxis().SetDecimals()
        graph_double_ratio.Draw("2pz")
        canvas_double_ratio.SaveAs(cfg["output"].replace(".root", "_double_ratio.pdf"))

    outfile = ROOT.TFile(cfg["output"], "recreate")
    graph_ratio_b.Write()
    graph_mid_b.Write()
    graph_fwd_b.Write()
    graph_mid_c.Write()
    graph_fwd_c.Write()
    graph_ratio_c.Write()
    if graph_double_ratio is not None:
        graph_double_ratio.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-c", metavar="text", default="config_fit.yml",
                        help="yaml config file for fit", required=True)
    args = parser.parse_args()
    plot(args.config)

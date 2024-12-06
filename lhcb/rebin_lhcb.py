"""
Script to rebin LHCb B+ cross section
"""

import sys
import argparse
import math
import numpy as np
import ROOT

sys.path.append('../utils') # pylint: disable=wrong-import-position
from style_formatter import set_global_style, set_object_style # pylint: disable=import-error


def check_bin_consistency(hist, ptlimits):
    """
    """

    ptlimits_orig = []
    for ipt in range(1, hist.GetNbinsX()+1):
        ptlimits_orig.append(hist.GetBinLowEdge(ipt))

    ptlimits_orig.append(hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX()))

    ptbins_ok = all(any(
        math.isclose(pt, pt_orig, rel_tol=1e-6) for pt_orig in ptlimits_orig) for pt in ptlimits)

    return ptbins_ok


# pylint: disable=too-many-arguments,too-many-locals
def get_rebinned_histos(hist, hist_stat, hist_syst, ptlimits, deltay, hist_names):
    """
    Function to rebin histogram for statistical uncertainties
    """

    if not check_bin_consistency(hist, ptlimits):
        print("ERROR: pt bins for rebinning not consistent with original histograms, exit")
        sys.exit()

    nb_to_pb = 1.e3

    hist_rebin_stat = ROOT.TH1F(
        hist_names["stat"],
        ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T}d#it{y} (pb #it{c}/GeV)",
        len(ptlimits)-1,
        np.array(ptlimits, np.float64)
    )
    hist_rebin_syst = ROOT.TH1F(
        hist_names["syst"],
        ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T}d#it{y} (pb #it{c}/GeV)",
        len(ptlimits)-1,
        np.array(ptlimits, np.float64)
    )
    set_object_style(hist_rebin_stat)
    set_object_style(hist_rebin_syst)
    for ipt, (pt_min, pt_max) in enumerate(zip(ptlimits[:-1], ptlimits[1:])):
        delta_pt = pt_max - pt_min
        cross_sec, cross_sec_stat, cross_sec_syst = 0., 0., 0.
        for ipt_orig in range(1, hist.GetNbinsX()+1):
            pt_low = hist.GetBinLowEdge(ipt_orig)
            pt_high = hist.GetXaxis().GetBinUpEdge(ipt_orig)
            if (pt_low > pt_min or math.isclose(pt_low, pt_min, rel_tol=1e-6)) and \
                (pt_high < pt_max or math.isclose(pt_high, pt_max, rel_tol=1e-6)):
                bin_width = hist.GetBinWidth(ipt_orig)
                cross_sec += hist.GetBinContent(ipt_orig) * bin_width
                cross_sec_stat += hist_stat.GetBinContent(ipt_orig)**2 * bin_width**2 # uncorrelated
                cross_sec_syst += hist_syst.GetBinContent(ipt_orig) * bin_width # correlated
        cross_sec = cross_sec / delta_pt / deltay * nb_to_pb
        cross_sec_stat = np.sqrt(cross_sec_stat) / delta_pt / deltay * nb_to_pb
        cross_sec_syst = cross_sec_syst / delta_pt / deltay * nb_to_pb
        hist_rebin_stat.SetBinContent(ipt+1, cross_sec)
        hist_rebin_syst.SetBinContent(ipt+1, cross_sec)
        hist_rebin_stat.SetBinError(ipt+1, cross_sec_stat)
        hist_rebin_syst.SetBinError(ipt+1, cross_sec_syst)

    return [hist_rebin_stat, hist_rebin_syst]

def rebin(infile_name, ptlimits, outfile_name):
    """
    Main function to rebin published LHCb results
    """

    set_global_style(padleftmargin=0.18, padbottommargin=0.14)

    infile = ROOT.TFile.Open(infile_name)
    inputdir = infile.Get("Table 3") # rapidity integrated measurement
    inputdir_rap = infile.Get("Table 1") # rapidity differential measurements
    hist_cent_rapint = inputdir.Get("Hist1D_y1")
    hist_stat_rapint = inputdir.Get("Hist1D_y1_e1")
    hist_syst_rapint = inputdir.Get("Hist1D_y1_e2")

    hist_reb_rapint = get_rebinned_histos(hist_cent_rapint, hist_stat_rapint,
                                          hist_syst_rapint, ptlimits, 2.5,
                                          {"stat": "h_crosssec_y_20_45_stat",
                                           "syst": "h_crosssec_y_20_45_syst"})

    hist_reb_rap = []
    for irap in range(5):
        hist_cent_rapint = inputdir_rap.Get(f"Hist1D_y{irap+1}")
        hist_stat_rapint = inputdir_rap.Get(f"Hist1D_y{irap+1}_e1")
        hist_syst_rapint = inputdir_rap.Get(f"Hist1D_y{irap+1}_e2")
        rap_min = 2.0 + irap * 0.5
        rap_max = 2.5 + irap * 0.5
        hist_reb_rap.append(get_rebinned_histos(hist_cent_rapint, hist_stat_rapint,
                                                hist_syst_rapint, ptlimits, 1., # already ds/dptdy
                                                {"stat": f"h_crosssec_y_{rap_min*10:.0f}_{rap_max*10:.0f}_stat",
                                                 "syst": f"h_crosssec_y_{rap_min*10:.0f}_{rap_max*10:.0f}_syst"}))

    outfile = ROOT.TFile(outfile_name, "recreate")
    hist_reb_rapint[0].Write()
    hist_reb_rapint[1].Write()
    for irap in range(5):
        hist_reb_rap[irap][0].Write()
        hist_reb_rap[irap][1].Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile", "-i", metavar="text",
                        default="HEPData-ins1630633-v1-root.root",
                        help="input HepData file for LHCb B+ measurement at 13 TeV")
    parser.add_argument("--ptlimits", "-p", nargs="+",
                        default=[2., 4., 6., 8., 10., 14., 23.5], type=float,
                        help="pt bins to rebin")
    parser.add_argument("--outfile", "-o", metavar="text",
                        default="LHCb_Bplus_13TeV_rebinned.root",
                        help="output file name")
    args = parser.parse_args()

    rebin(args.infile, args.ptlimits, args.outfile)

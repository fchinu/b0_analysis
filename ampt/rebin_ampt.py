import math
import sys
import numpy as np
import ROOT
sys.path.append('utils') # pylint: disable=wrong-import-position
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

def get_rebinned_histos(hist, ptlimits, hist_name):
    """
    Function to rebin histogram for statistical uncertainties
    """

    if not check_bin_consistency(hist, ptlimits):
        print("ERROR: pt bins for rebinning not consistent with original histograms, exit")
        sys.exit()

    hist_rebin_stat = ROOT.TH1F(
        hist_name,
        ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T}d#it{y} (#mub #it{c}/GeV)",
        len(ptlimits)-1,
        np.array(ptlimits, np.float64)
    )
    set_object_style(hist_rebin_stat)
    for ipt, (pt_min, pt_max) in enumerate(zip(ptlimits[:-1], ptlimits[1:])):
        delta_pt = pt_max - pt_min
        cross_sec, cross_sec_stat = 0., 0.
        for ipt_orig in range(1, hist.GetNbinsX()+1):
            pt_low = hist.GetBinLowEdge(ipt_orig)
            pt_high = hist.GetXaxis().GetBinUpEdge(ipt_orig)
            if (pt_low > pt_min or math.isclose(pt_low, pt_min, rel_tol=1e-6)) and \
                (pt_high < pt_max or math.isclose(pt_high, pt_max, rel_tol=1e-6)):
                bin_width = hist.GetBinWidth(ipt_orig)
                cross_sec += hist.GetBinContent(ipt_orig) * bin_width
                cross_sec_stat += hist.GetBinError(ipt_orig)**2 * bin_width**2 # uncorrelated
        cross_sec = cross_sec / delta_pt
        cross_sec_stat = np.sqrt(cross_sec_stat) / delta_pt
        hist_rebin_stat.SetBinContent(ipt+1, cross_sec)
        hist_rebin_stat.SetBinError(ipt+1, cross_sec_stat)

    return hist_rebin_stat

if __name__ == "__main__":
    ampt_file_48 = ROOT.TFile.Open("ampt/AMPT_13TeV_mb4.8.root")
    ampt_file_66 = ROOT.TFile.Open("ampt/AMPT_13TeV_mb6.6.root")

    histos_to_rebin = ["B0_mid", "D0_mid",
                       "Bp_forward_2_25", "Bp_forward_3_35", "Bp_forward_4_45",
                       "B0_forward_2_25", "B0_forward_3_35", "B0_forward_4_45",
                       "D0_forward_2_25", "D0_forward_3_35", "D0_forward_4_45",
                      ]
    
    rebinned_histos_48 = {}
    rebinned_histos_66 = {}
    set_global_style()

    for histo_name in histos_to_rebin:
        for file in [ampt_file_48, ampt_file_66]:
            hist = file.Get(histo_name)
            pt_bins = [0., 1., 2., 4., 6., 8., 10., 14.]
            if file == ampt_file_48:
                rebinned_histos_48[histo_name] = get_rebinned_histos(hist, pt_bins, histo_name+"_rebin")
            else:
                rebinned_histos_66[histo_name] = get_rebinned_histos(hist, pt_bins, histo_name+"_rebin")

    ratio_histos_names = [
        "B0_over_Bp_2_25", "B0_over_Bp_3_35", "B0_over_Bp_4_45",
        "D0_mid_over_forward_2_25", "D0_mid_over_forward_3_35", "D0_mid_over_forward_4_45",
    ]
    ratio_histos_48 = {}
    ratio_histos_66 = {}

    for ratio_name in ratio_histos_names:
        if "B0_over_Bp" in ratio_name:
            rap = ratio_name.split("B0_over_Bp_")[-1]
            hist_b0_48 = rebinned_histos_48["B0_mid"]
            hist_bp_48 = rebinned_histos_48[f"Bp_forward_{rap}"]
            hist_b0_66 = rebinned_histos_66["B0_mid"]
            hist_bp_66 = rebinned_histos_66[f"Bp_forward_{rap}"]
            ratio_histos_48[ratio_name] = hist_b0_48.Clone(ratio_name+"_13TeV_mb4.8")
            ratio_histos_48[ratio_name].Divide(hist_bp_48)
            ratio_histos_66[ratio_name] = hist_b0_66.Clone(ratio_name+"_13TeV_mb6.6")
            ratio_histos_66[ratio_name].Divide(hist_bp_66)
        elif "D0_mid_over_forward" in ratio_name:
            rap = ratio_name.split("D0_mid_over_forward_")[-1]
            hist_d0_mid_48 = rebinned_histos_48["D0_mid"]
            hist_d0_fwd_48 = rebinned_histos_48[f"D0_forward_{rap}"]
            hist_d0_mid_66 = rebinned_histos_66["D0_mid"]
            hist_d0_fwd_66 = rebinned_histos_66[f"D0_forward_{rap}"]
            ratio_histos_48[ratio_name] = hist_d0_mid_48.Clone(ratio_name+"_13TeV_mb4.8")
            ratio_histos_48[ratio_name].Divide(hist_d0_fwd_48)
            ratio_histos_66[ratio_name] = hist_d0_mid_66.Clone(ratio_name+"_13TeV_mb6.6")
            ratio_histos_66[ratio_name].Divide(hist_d0_fwd_66)

    outfile_48 = ROOT.TFile("ampt/AMPT_13TeV_mb4.8_rebinned.root", "recreate")
    for histo in rebinned_histos_48.values():
        histo.Write()
    for histo in ratio_histos_48.values():
        histo.Write()
    outfile_48.Close()
    outfile_66 = ROOT.TFile("ampt/AMPT_13TeV_mb6.6_rebinned.root", "recreate")
    for histo in rebinned_histos_66.values():
        histo.Write()
    for histo in ratio_histos_66.values():
        histo.Write()
    outfile_66.Close()
    ampt_file_48.Close()
    ampt_file_66.Close()

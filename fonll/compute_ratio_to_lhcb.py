"""
"""

import argparse
import sys
sys.path.append('../utils') # pylint: disable=wrong-import-position
from style_formatter import set_global_style, set_object_style # pylint: disable=import-error
import ROOT


def compute_ratio(infilename_alice, infilename_lhcb): # pylint: disable=too-many-locals
    """
    Main function for ratio calculation
    """

    set_global_style(padbottommargin=0.14, padleftmargin=0.16)

    infile_alice = ROOT.TFile.Open(infilename_alice)
    hist_alice_stat = infile_alice.Get("h_stat")
    hist_alice_syst = infile_alice.Get("h_syst")
    set_object_style(hist_alice_stat, markersize=1.5)
    set_object_style(hist_alice_syst, fillstyle=0)

    infile_lhcb = ROOT.TFile.Open(infilename_lhcb)
    rapidities_lhcb  = ["20_25", "30_35", "40_45"]
    hist_lhcb_stat, hist_lhcb_syst, hist_ratio_stat, hist_ratio_syst = ([] for _ in range(4))
    for irap, rapidity in enumerate(rapidities_lhcb):
        hist_lhcb_stat.append(infile_lhcb.Get(f"h_crosssec_y_{rapidity}_stat"))
        hist_lhcb_syst.append(infile_lhcb.Get(f"h_crosssec_y_{rapidity}_syst"))
        set_object_style(hist_lhcb_stat[-1])
        set_object_style(hist_lhcb_syst[-1], fillstyle=0)
        hist_ratio_stat.append(hist_alice_stat.Clone(f"hist_ratio_stat_{rapidity}"))
        hist_ratio_stat[-1].Divide(hist_lhcb_stat[-1])
        hist_ratio_syst.append(hist_alice_syst.Clone(f"hist_ratio_syst_{rapidity}"))
        hist_ratio_syst[-1].Divide(hist_lhcb_syst[-1])

    graph_ratio_fonll = []
    for irap, rapidity in enumerate(rapidities_lhcb):
        infile_fonll = ROOT.TFile.Open(f"fonll_bhadron_nnpdfs_y05_y{rapidity}_13tev_13dot6tev.root")
        graph_ratio_fonll.append(infile_fonll.Get("graph_ratio_mid_fwd"))
        set_object_style(graph_ratio_fonll[-1], color=ROOT.kAzure+4, alpha=0.5)

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.045)
    lat.SetTextFont(42)
    lat.SetTextColor(ROOT.kBlack)

    leg_ratio = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
    leg_ratio.SetBorderSize(0)
    leg_ratio.SetFillStyle(0)
    leg_ratio.SetTextSize(0.045)
    leg_ratio.AddEntry(hist_ratio_syst[-1], "Data", "lp")
    leg_ratio.AddEntry(graph_ratio_fonll[-1], "FONLL", "f")

    line_at_one = ROOT.TLine(0., 1., 28., 1.)
    line_at_one.SetLineWidth(1)
    line_at_one.SetLineStyle(9)
    line_at_one.SetLineColor(ROOT.kGray+2)

    labels_lhcb = ["LHCb #sqrt{#it{s}} = 13 TeV, 2.0 < #it{y} < 2.5",
                   "LHCb #sqrt{#it{s}} = 13 TeV, 3.0 < #it{y} < 3.5",
                   "LHCb #sqrt{#it{s}} = 13 TeV, 4.0 < #it{y} < 4.5"]

    canv_ratio = ROOT.TCanvas("canv_ratio", "", 1500, 500)
    canv_ratio.Divide(3, 1)
    for irap, label_lhcb in enumerate(labels_lhcb):
        canv_ratio.cd(irap + 1).SetLogy()
        hframe = canv_ratio.cd(irap + 1).DrawFrame(
            0., 0.3, 28., 20.,
            ";#it{p}_{T} (GeV/#it{c}); B^{0} (mid) / B^{#plus} (fwd)"
        )
        hframe.GetYaxis().SetMoreLogLabels()
        line_at_one.Draw()
        graph_ratio_fonll[irap].Draw("2")
        hist_ratio_syst[irap].DrawCopy("e2same")
        hist_ratio_stat[irap].DrawCopy("same")
        lat.DrawLatex(0.2, 0.26, "ALICE #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5")
        lat.DrawLatex(0.2, 0.2, label_lhcb)
        if irap == 0:
            leg_ratio.Draw()
    canv_ratio.SaveAs("B_mid_fwd_ratio.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile_alice", "-a", metavar="text",
                        default="../cross_section_default_w_syst.root",
                        help="input file with ALICE results", required=False)
    parser.add_argument("--infile_lhcb", "-b", metavar="text",
                        default="../lhcb/LHCb_Bplus_13TeV_rebinned.root",
                        help="input file with LHCb results", required=False)
    args = parser.parse_args()

    compute_ratio(args.infile_alice, args.infile_lhcb)

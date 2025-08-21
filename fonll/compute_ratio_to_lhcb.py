"""
"""

import argparse
import sys
sys.path.append('../utils') # pylint: disable=wrong-import-position
from style_formatter import set_global_style, set_object_style # pylint: disable=import-error
import ROOT

def convert_to_graph(histo):
    """
    Function to convert histo to graph
    """

    nbins = histo.GetNbinsX()
    graph = ROOT.TGraphErrors(nbins)
    graph.SetName(f"graph_from_{histo.GetName()}")
    set_object_style(graph, fillstyle=0, color=histo.GetMarkerColor())
    for ipt in range(nbins):
        graph.SetPoint(ipt, histo.GetBinCenter(ipt+1), histo.GetBinContent(ipt+1))
        graph.SetPointError(ipt, histo.GetBinWidth(ipt+1)/4, histo.GetBinError(ipt+1))

    return graph


def compute_ratio(infilename_alice, infilename_lhcb): # pylint: disable=too-many-locals
    """
    Main function for ratio calculation
    """

    set_global_style(padbottommargin=0.14, padleftmargin=0.16, titlesize=0.05, lablesize=0.05)

    # ALICE B0
    infile_alice = ROOT.TFile.Open(infilename_alice)
    hist_alice_stat = infile_alice.Get("h_stat")
    hist_alice_syst = infile_alice.Get("h_syst")
    set_object_style(hist_alice_stat, markersize=1.5)
    set_object_style(hist_alice_syst, fillstyle=0)

    # LHCb B+
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

    # ALICE/LHCb D0
    hist_ratio_d_stat, graph_ratio_d_syst = [], []
    infile_alice_d = ROOT.TFile.Open("HEPData-ins2697877-v1-root.root")
    for irap in range(3):
        hist_stat_unc = infile_alice_d.Get(f"Figure 10 D0 y range {irap}/Hist1D_y1_e1")
        hist_syst_unc_low = infile_alice_d.Get(f"Figure 10 D0 y range {irap}/Hist1D_y1_e2minus")
        hist_syst_unc_high = infile_alice_d.Get(f"Figure 10 D0 y range {irap}/Hist1D_y1_e2plus")
        hist_ratio_d_stat.append(infile_alice_d.Get(f"Figure 10 D0 y range {irap}/Hist1D_y1"))
        graph_ratio_d_syst.append(ROOT.TGraphAsymmErrors(hist_stat_unc.GetNbinsX()))
        for ipt in range(1, hist_ratio_d_stat[irap].GetNbinsX()+1):
            hist_ratio_d_stat[irap].SetBinError(ipt, hist_stat_unc.GetBinContent(ipt))
            graph_ratio_d_syst[irap].SetPoint(ipt-1, hist_ratio_d_stat[irap].GetBinCenter(ipt),
                                              hist_ratio_d_stat[irap].GetBinContent(ipt))
            pt_unc = hist_ratio_d_stat[irap].GetBinWidth(ipt)/4
            graph_ratio_d_syst[irap].SetPointError(ipt-1, pt_unc, pt_unc,
                                                   -hist_syst_unc_low.GetBinContent(ipt),
                                                   hist_syst_unc_high.GetBinContent(ipt))
        set_object_style(hist_ratio_d_stat[irap], color=ROOT.kBlack)
        set_object_style(graph_ratio_d_syst[irap], color=ROOT.kBlack, fillstyle=0)

    # FONLL
    graph_ratio_fonll, graph_ratiotocent_fonll = [], []
    for irap, rapidity in enumerate(rapidities_lhcb):
        infile_fonll = ROOT.TFile.Open(
            f"fonll_bhadron_chadron_nnpdfs_y05_y{rapidity}_13tev.root")
        graph_ratio_fonll.append(infile_fonll.Get("graph_ratio_mid_fwd_b"))
        set_object_style(graph_ratio_fonll[irap], color=ROOT.kAzure+4, alpha=0.5)
        graph_ratiotocent_fonll.append(graph_ratio_fonll[-1].Clone())
        for ipt in range(graph_ratio_fonll[irap].GetN()):
            ratio = graph_ratio_fonll[irap].GetPointY(ipt)
            ratio_unchigh = ratio + graph_ratio_fonll[irap].GetErrorYhigh(ipt)
            ratio_unclow = ratio - graph_ratio_fonll[irap].GetErrorYlow(ipt)
            graph_ratiotocent_fonll[irap].SetPoint(ipt, graph_ratio_fonll[-1].GetPointX(ipt), 1.)
            graph_ratiotocent_fonll[irap].SetPointError(
                ipt, graph_ratio_fonll[irap].GetErrorXlow(ipt),
                graph_ratio_fonll[irap].GetErrorXhigh(ipt),
                1 - ratio_unclow / ratio,
                ratio_unchigh/ratio - 1)

    graph_ratio_fonll_onlycent = []
    for irap, _ in enumerate(rapidities_lhcb):
        graph_ratio_fonll_onlycent.append(graph_ratio_fonll[irap].Clone(
            f"graph_ratio_fonll_onlycent_{irap}"))
        for ipt in range(graph_ratio_fonll_onlycent[irap].GetN()):
            graph_ratio_fonll_onlycent[irap].SetPointError(
                ipt,
                graph_ratio_fonll[irap].GetErrorXlow(ipt),
                graph_ratio_fonll[irap].GetErrorXhigh(ipt),
                0., 0.)

            set_object_style(graph_ratio_fonll_onlycent[-1], color=ROOT.kAzure+4, markersize=0)

    graph_ratio_d_fonll, graph_ratiotocent_d_fonll = [], []
    for irap, rapidity in enumerate(rapidities_lhcb):
        infile_fonll = ROOT.TFile.Open(
            f"fonll_bhadron_chadron_nnpdfs_y05_y{rapidity}_13tev.root")
        graph_ratio_d_fonll.append(infile_fonll.Get("graph_ratio_mid_fwd_c"))
        set_object_style(graph_ratio_d_fonll[irap], color=ROOT.kRed+1, alpha=0.5)
        graph_ratiotocent_d_fonll.append(graph_ratio_d_fonll[-1].Clone())
        for ipt in range(graph_ratio_d_fonll[irap].GetN()):
            ratio = graph_ratio_d_fonll[irap].GetPointY(ipt)
            ratio_unchigh = ratio + graph_ratio_d_fonll[irap].GetErrorYhigh(ipt)
            ratio_unclow = ratio - graph_ratio_d_fonll[irap].GetErrorYlow(ipt)
            graph_ratiotocent_d_fonll[irap].SetPoint(ipt, graph_ratio_d_fonll[-1].GetPointX(ipt), 1.)
            graph_ratiotocent_d_fonll[irap].SetPointError(
                ipt, graph_ratio_d_fonll[irap].GetErrorXlow(ipt),
                graph_ratio_d_fonll[irap].GetErrorXhigh(ipt),
                1 - ratio_unclow / ratio,
                ratio_unchigh/ratio - 1)

    graph_ratio_d_fonll_onlycent = []
    for irap, _ in enumerate(rapidities_lhcb):
        graph_ratio_d_fonll_onlycent.append(graph_ratio_d_fonll[irap].Clone(
            f"graph_ratio_d_fonll_onlycent_{irap}"))
        for ipt in range(graph_ratio_d_fonll_onlycent[irap].GetN()):
            graph_ratio_d_fonll_onlycent[irap].SetPointError(
                ipt,
                graph_ratio_d_fonll[irap].GetErrorXlow(ipt),
                graph_ratio_d_fonll[irap].GetErrorXhigh(ipt),
                0., 0.)

            set_object_style(graph_ratio_d_fonll_onlycent[-1], color=ROOT.kRed+1, markersize=0)

    # Data / FONLL
    hist_ratiotofonll_stat, hist_ratiotofonll_syst = [], []
    for irap in range(3):
        hist_ratiotofonll_stat.append(hist_ratio_stat[irap].Clone())
        hist_ratiotofonll_syst.append(hist_ratio_syst[irap].Clone())
        for ipt in range(hist_ratiotofonll_stat[irap].GetNbinsX()):
            ratio_fonll = graph_ratio_fonll[irap].GetPointY(ipt)
            ratio_data = hist_ratiotofonll_stat[irap].GetBinContent(ipt+1)
            ratio_data_stat = hist_ratiotofonll_stat[irap].GetBinError(ipt+1)
            ratio_data_syst = hist_ratiotofonll_syst[irap].GetBinError(ipt+1)
            hist_ratiotofonll_stat[irap].SetBinContent(ipt+1, ratio_data / ratio_fonll)
            hist_ratiotofonll_stat[irap].SetBinError(ipt+1, ratio_data_stat / ratio_fonll)
            hist_ratiotofonll_syst[irap].SetBinContent(ipt+1, ratio_data / ratio_fonll)
            hist_ratiotofonll_syst[irap].SetBinError(ipt+1, ratio_data_syst / ratio_fonll)

    # Data / FONLL
    hist_ratiotofonll_d_stat, graph_ratiotofonll_d_syst = [], []
    for irap in range(3):
        hist_ratiotofonll_d_stat.append(hist_ratio_d_stat[irap].Clone())
        graph_ratiotofonll_d_syst.append(graph_ratio_d_syst[irap].Clone())
        for ipt in range(hist_ratiotofonll_d_stat[irap].GetNbinsX()):
            ratio_fonll = graph_ratio_d_fonll[irap].GetPointY(ipt)
            ratio_data = hist_ratiotofonll_d_stat[irap].GetBinContent(ipt+1)
            ratio_data_stat = hist_ratiotofonll_d_stat[irap].GetBinError(ipt+1)
            ratio_data_systhigh = graph_ratiotofonll_d_syst[irap].GetErrorYhigh(ipt)
            ratio_data_systlow = graph_ratiotofonll_d_syst[irap].GetErrorYlow(ipt)
            hist_ratiotofonll_d_stat[irap].SetBinContent(ipt+1, ratio_data / ratio_fonll)
            hist_ratiotofonll_d_stat[irap].SetBinError(ipt+1, ratio_data_stat / ratio_fonll)
            graph_ratiotofonll_d_syst[irap].SetPoint(ipt, hist_ratiotofonll_d_stat[irap].GetBinCenter(ipt+1),
                                                     ratio_data / ratio_fonll)
            graph_ratiotofonll_d_syst[irap].SetPointError(ipt, graph_ratiotofonll_d_syst[irap].GetErrorXlow(ipt),
                                                          graph_ratiotofonll_d_syst[irap].GetErrorXhigh(ipt),
                                                          ratio_data_systlow / ratio_fonll,
                                                          ratio_data_systhigh / ratio_fonll)

    lat = [ROOT.TLatex(), ROOT.TLatex(), ROOT.TLatex()]
    size_factor = [1., 0.36/0.315, 0.36/0.325]
    for irap in range(3):
        lat[irap].SetNDC()
        lat[irap].SetTextSize(0.05 * size_factor[irap])
        lat[irap].SetTextFont(42)
        lat[irap].SetTextColor(ROOT.kBlack)

    lat_alice = ROOT.TLatex()
    lat_alice.SetNDC()
    lat_alice.SetTextSize(0.065)
    lat_alice.SetTextFont(42)
    lat_alice.SetTextColor(ROOT.kBlack)

    lat_labels = ROOT.TLatex()
    lat_labels.SetNDC()
    lat_labels.SetTextSize(0.055)
    lat_labels.SetTextFont(42)
    lat_labels.SetTextColor(ROOT.kBlack)

    leg_ratio = ROOT.TLegend(0.04, 0.7, 0.5, 0.9)
    leg_ratio.SetBorderSize(0)
    leg_ratio.SetFillStyle(0)
    leg_ratio.SetTextSize(0.05 * size_factor[1])
    leg_ratio.AddEntry(hist_ratio_syst[-1], "#frac{B^{0} ALICE (13.6 TeV)}{B^{#plus} LHCb (13 TeV)}", "lp")
    leg_ratio.AddEntry(graph_ratio_fonll[-1], "FONLL", "f")

    leg_ratio_cent = ROOT.TLegend(0.04, 0.7, 0.5, 0.9)
    leg_ratio_cent.SetBorderSize(0)
    leg_ratio_cent.SetFillStyle(0)
    leg_ratio_cent.SetTextSize(0.05 * size_factor[1])
    leg_ratio_cent.AddEntry(hist_ratio_syst[-1], "#frac{B^{0} ALICE (13.6 TeV)}{B^{#plus} LHCb (13 TeV)}", "lp")
    leg_ratio_cent.AddEntry(graph_ratio_fonll_onlycent[-1], "FONLL", "l")

    leg_ratio_d = ROOT.TLegend(0.08, 0.7, 0.54, 0.9)
    leg_ratio_d.SetBorderSize(0)
    leg_ratio_d.SetFillStyle(0)
    leg_ratio_d.SetTextSize(0.05 * size_factor[1])
    leg_ratio_d.AddEntry(hist_ratio_d_stat[-1], "#frac{D^{0} ALICE (13 TeV)}{D^{0} LHCb (13 TeV)}", "lp")
    leg_ratio_d.AddEntry(graph_ratio_d_fonll[-1], "FONLL", "f")

    leg_ratio_d_cent = ROOT.TLegend(0.08, 0.7, 0.54, 0.9)
    leg_ratio_d_cent.SetBorderSize(0)
    leg_ratio_d_cent.SetFillStyle(0)
    leg_ratio_d_cent.SetTextSize(0.05 * size_factor[1])
    leg_ratio_d_cent.AddEntry(hist_ratio_d_stat[-1], "#frac{D^{0} ALICE (13 TeV)}{D^{0} LHCb (13 TeV)}", "lp")
    leg_ratio_d_cent.AddEntry(graph_ratio_d_fonll_onlycent[-1], "FONLL", "l")

    line_at_one = ROOT.TLine(0., 1., 27., 1.)
    line_at_one.SetLineWidth(1)
    line_at_one.SetLineStyle(9)
    line_at_one.SetLineColor(ROOT.kGray+2)

    labels_lhcb = ["LHCb 2.0 < #it{y} < 2.5",
                   "LHCb 3.0 < #it{y} < 3.5",
                   "LHCb 4.0 < #it{y} < 4.5"]

    canv_ratio = ROOT.TCanvas("canv_ratio", "", 1400, 600)
    pad_ratio = [ROOT.TPad("pad_ratio_1", "pad_ratio_1", 0., 0., 0.36, 1.),
                 ROOT.TPad("pad_ratio_2", "pad_ratio_2", 0.36, 0., 0.675, 1.),
                 ROOT.TPad("pad_ratio_3", "pad_ratio_3", 0.675, 0., 1., 1.)]
    pad_ratio[0].SetRightMargin(0)
    pad_ratio[1].SetLeftMargin(0)
    pad_ratio[1].SetRightMargin(0)
    pad_ratio[2].SetLeftMargin(0)

    graph_ratio_syst = [None, None, None]
    for irap in [2, 1, 0]:
        canv_ratio.cd()
        pad_ratio[irap].Draw()
        hframe = pad_ratio[irap].cd().DrawFrame(
            0.01, 0.3, 27., 20.,
            ";#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} mid / forward rapidity"
        )
        pad_ratio[irap].cd().SetLogy()
        hframe.GetYaxis().SetMoreLogLabels()
        hframe.GetYaxis().SetNoExponent()
        hframe.GetYaxis().SetTitleSize(0.05 * size_factor[irap])
        hframe.GetXaxis().SetTitleSize(0.05 * size_factor[irap])
        hframe.GetYaxis().SetLabelSize(0.05 * size_factor[irap])
        hframe.GetXaxis().SetLabelSize(0.05 * size_factor[irap])
        hframe.GetXaxis().SetTitleOffset(1.13 * 1./size_factor[irap])
        line_at_one.Draw()
        graph_ratio_fonll[irap].Draw("2")
        graph_ratio_fonll_onlycent[irap].Draw("pz")
        # graph_ratio_d_syst[irap].Draw("2")
        # hist_ratio_d_stat[irap].DrawCopy("same")
        graph_ratio_syst[irap] = convert_to_graph(hist_ratio_syst[irap])
        graph_ratio_syst[irap].Draw("2")
        hist_ratio_stat[irap].DrawCopy("same")
        xmin = 0.2
        if irap > 0:
            xmin = 0.04
        else:
            lat_alice.DrawLatex(xmin+0.02, 0.88, "ALICE")
            lat_labels.DrawLatex(xmin+0.02, 0.82, "pp collisions")
        lat[irap].DrawLatex(xmin, 0.26, "ALICE |#it{y}| < 0.5")
        lat[irap].DrawLatex(xmin, 0.2, labels_lhcb[irap])
        if irap == 1:
            leg_ratio.Draw()
            leg_ratio_cent.Draw()
    canv_ratio.SaveAs("B_mid_fwd_ratio_vsFONLL.pdf")

    canv_ratio_lin = ROOT.TCanvas("canv_ratio_lin", "", 1400, 600)
    canv_ratio_lin.Divide(3, 1)
    max_ratio = [2., 4., 15.]

    for irap in [2, 1, 0]:
        hframe = canv_ratio_lin.cd(irap+1).DrawFrame(
            0.01, 0., 27., max_ratio[irap],
            ";#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} mid / forward rapidity"
        )
        hframe.GetYaxis().SetDecimals()
        line_at_one.Draw()
        graph_ratio_fonll[irap].Draw("2")
        graph_ratio_fonll_onlycent[irap].Draw("pz")
        # graph_ratio_d_syst[irap].Draw("2")
        # hist_ratio_d_stat[irap].DrawCopy("same")
        graph_ratio_syst[irap].Draw("2")
        hist_ratio_stat[irap].DrawCopy("same")
        shift_x = 0
        if irap == 2:
            shift_x += 0.25
        lat[0].DrawLatex(0.2 + shift_x, 0.26, "ALICE |#it{y}| < 0.5")
        lat[0].DrawLatex(0.2 + shift_x, 0.22, labels_lhcb[irap])
        if irap == 1:
            leg_ratio.SetX1NDC(0.2)
            leg_ratio.SetX2NDC(0.7)
            leg_ratio_cent.SetX1NDC(0.2)
            leg_ratio_cent.SetX2NDC(0.7)
            leg_ratio.Draw()
            leg_ratio_cent.Draw()
    canv_ratio_lin.SaveAs("B_mid_fwd_ratio_vsFONLL_linscale.pdf")

    canv_ratio_wratio = ROOT.TCanvas("canv_ratio_wratio", "", 1400, 1200)
    pad_ratio_wratio = [ROOT.TPad("pad_ratio_wratio_1", "pad_ratio_wratio_1", 0., 0.44, 0.36, 1.0),
                        ROOT.TPad("pad_ratio_wratio_2", "pad_ratio_wratio_2", 0.36, 0.44, 0.675, 1.0),
                        ROOT.TPad("pad_ratio_wratio_3", "pad_ratio_wratio_3", 0.675, 0.44, 1., 1.0),
                        ROOT.TPad("pad_ratio_wratio_4", "pad_ratio_wratio_4", 0., 0., 0.36, 0.52),
                        ROOT.TPad("pad_ratio_wratio_5", "pad_ratio_wratio_5", 0.36, 0., 0.675, 0.52),
                        ROOT.TPad("pad_ratio_wratio_6", "pad_ratio_wratio_6", 0.675, 0., 1., 0.52)]
    pad_ratio_wratio[0].SetRightMargin(0)
    pad_ratio_wratio[1].SetLeftMargin(0)
    pad_ratio_wratio[1].SetRightMargin(0)
    pad_ratio_wratio[2].SetLeftMargin(0)
    pad_ratio_wratio[3].SetRightMargin(0)
    pad_ratio_wratio[4].SetLeftMargin(0)
    pad_ratio_wratio[4].SetRightMargin(0)
    pad_ratio_wratio[5].SetLeftMargin(0)
    pad_ratio_wratio[3].SetTopMargin(0)
    pad_ratio_wratio[4].SetTopMargin(0)
    pad_ratio_wratio[5].SetTopMargin(0)

    graph_ratiotofonll_syst = [None, None, None, None, None, None]
    for irap in [2, 1, 0, 5, 4, 3]:
        canv_ratio_wratio.cd()
        pad_ratio_wratio[irap].Draw()
        if irap < 3:
            hframe = pad_ratio_wratio[irap].cd().DrawFrame(
                0.01, 0.24, 27., 20.,
                ";#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} mid / forward rapidity"
            )
            pad_ratio_wratio[irap].cd().SetLogy()
            hframe.GetYaxis().SetMoreLogLabels()
            hframe.GetYaxis().SetNoExponent()
            hframe.GetYaxis().SetTitleSize(0.05 * size_factor[irap])
            hframe.GetXaxis().SetTitleSize(0.05 * size_factor[irap])
            hframe.GetYaxis().SetLabelSize(0.05 * size_factor[irap])
            hframe.GetXaxis().SetLabelSize(0.05 * size_factor[irap])
            hframe.GetXaxis().SetTitleOffset(1.13 * 1./size_factor[irap])
            line_at_one.Draw()
            graph_ratio_fonll[irap].Draw("2")
            graph_ratio_fonll_onlycent[irap].Draw("pz")
            graph_ratio_syst[irap] = convert_to_graph(hist_ratio_syst[irap])
            graph_ratio_syst[irap].Draw("2")
            hist_ratio_stat[irap].DrawCopy("same")
            xmin = 0.2
            if irap > 0:
                xmin = 0.04
            else:
                lat_alice.DrawLatex(xmin+0.02, 0.88, "ALICE")
                lat_labels.DrawLatex(xmin+0.02, 0.82, "pp collisions")
            lat[irap].DrawLatex(xmin, 0.26, "ALICE |#it{y}| < 0.5")
            lat[irap].DrawLatex(xmin, 0.2, labels_lhcb[irap])
            if irap == 1:
                leg_ratio.Draw()
                leg_ratio_cent.Draw()
        else:
            hframe = pad_ratio_wratio[irap].cd().DrawFrame(
                0.01, 0., 27., 2.1,
                ";#it{p}_{T} (GeV/#it{c}); data / FONLL"
            )
            hframe.GetYaxis().SetDecimals()
            hframe.GetYaxis().SetTitleSize(0.05 * size_factor[irap-3])
            hframe.GetXaxis().SetTitleSize(0.05 * size_factor[irap-3])
            hframe.GetYaxis().SetLabelSize(0.05 * size_factor[irap-3])
            hframe.GetXaxis().SetLabelSize(0.05 * size_factor[irap-3])
            hframe.GetXaxis().SetTitleOffset(1.13 * 1./size_factor[irap-3])
            line_at_one.Draw()
            graph_ratiotocent_fonll[irap-3].Draw("2")
            graph_ratiotofonll_syst[irap-3] = convert_to_graph(hist_ratiotofonll_syst[irap-3])
            graph_ratiotofonll_syst[irap-3].Draw("2")
            hist_ratiotofonll_stat[irap-3].DrawCopy("same")
            hist_ratiotofonll_stat[irap-3].Fit("pol1", "0")

    canv_ratio_wratio.SaveAs("B_mid_fwd_ratio_vsFONLL_withRatioToFONLL.pdf")


    canv_ratio_d_wratio = ROOT.TCanvas("canv_ratio_d_wratio", "", 1400, 1200)
    pad_ratio_d_wratio = [ROOT.TPad("pad_ratio_d_wratio_1", "pad_ratio_d_wratio_1", 0., 0.44, 0.36, 1.0),
                          ROOT.TPad("pad_ratio_d_wratio_2", "pad_ratio_d_wratio_2", 0.36, 0.44, 0.675, 1.0),
                          ROOT.TPad("pad_ratio_d_wratio_3", "pad_ratio_d_wratio_3", 0.675, 0.44, 1., 1.0),
                          ROOT.TPad("pad_ratio_d_wratio_4", "pad_ratio_d_wratio_4", 0., 0., 0.36, 0.52),
                          ROOT.TPad("pad_ratio_d_wratio_5", "pad_ratio_d_wratio_5", 0.36, 0., 0.675, 0.52),
                          ROOT.TPad("pad_ratio_d_wratio_6", "pad_ratio_d_wratio_6", 0.675, 0., 1., 0.52)]
    pad_ratio_d_wratio[0].SetRightMargin(0)
    pad_ratio_d_wratio[1].SetLeftMargin(0)
    pad_ratio_d_wratio[1].SetRightMargin(0)
    pad_ratio_d_wratio[2].SetLeftMargin(0)
    pad_ratio_d_wratio[3].SetRightMargin(0)
    pad_ratio_d_wratio[4].SetLeftMargin(0)
    pad_ratio_d_wratio[4].SetRightMargin(0)
    pad_ratio_d_wratio[5].SetLeftMargin(0)
    pad_ratio_d_wratio[3].SetTopMargin(0)
    pad_ratio_d_wratio[4].SetTopMargin(0)
    pad_ratio_d_wratio[5].SetTopMargin(0)

    graph_ratiotofonll_syst = [None, None, None, None, None, None]
    for irap in [2, 1, 0, 5, 4, 3]:
        canv_ratio_d_wratio.cd()
        pad_ratio_d_wratio[irap].Draw()
        if irap < 3:
            hframe = pad_ratio_d_wratio[irap].cd().DrawFrame(
                0.01, 0.24, 14., 20.,
                ";#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} mid / forward rapidity"
            )
            pad_ratio_d_wratio[irap].cd().SetLogy()
            hframe.GetYaxis().SetMoreLogLabels()
            hframe.GetYaxis().SetNoExponent()
            hframe.GetYaxis().SetTitleSize(0.05 * size_factor[irap])
            hframe.GetXaxis().SetTitleSize(0.05 * size_factor[irap])
            hframe.GetYaxis().SetLabelSize(0.05 * size_factor[irap])
            hframe.GetXaxis().SetLabelSize(0.05 * size_factor[irap])
            hframe.GetXaxis().SetTitleOffset(1.13 * 1./size_factor[irap])
            line_at_one.Draw()
            graph_ratio_d_fonll[irap].Draw("2")
            graph_ratio_d_fonll_onlycent[irap].Draw("pz")
            graph_ratio_d_syst[irap].Draw("2")
            hist_ratio_d_stat[irap].DrawCopy("same")
            xmin = 0.3
            if irap > 0:
                xmin = 0.14
            else:
                lat_alice.DrawLatex(0.22, 0.88, "ALICE")
                lat_labels.DrawLatex(0.22, 0.82, "pp collisions")
            lat[irap].DrawLatex(xmin, 0.26, "ALICE |#it{y}| < 0.5")
            lat[irap].DrawLatex(xmin, 0.2, labels_lhcb[irap])
            if irap == 1:
                leg_ratio_d.Draw()
                leg_ratio_d_cent.Draw()
        else:
            hframe = pad_ratio_d_wratio[irap].cd().DrawFrame(
                0.01, 0., 14., 2.1,
                ";#it{p}_{T} (GeV/#it{c}); data / FONLL"
            )
            hframe.GetYaxis().SetDecimals()
            hframe.GetYaxis().SetTitleSize(0.05 * size_factor[irap-3])
            hframe.GetXaxis().SetTitleSize(0.05 * size_factor[irap-3])
            hframe.GetYaxis().SetLabelSize(0.05 * size_factor[irap-3])
            hframe.GetXaxis().SetLabelSize(0.05 * size_factor[irap-3])
            hframe.GetXaxis().SetTitleOffset(1.13 * 1./size_factor[irap-3])
            line_at_one.Draw()
            graph_ratiotocent_d_fonll[irap-3].Draw("2")
            graph_ratiotofonll_d_syst[irap-3].Draw("2")
            hist_ratiotofonll_d_stat[irap-3].DrawCopy("same")
            hist_ratiotofonll_d_stat[irap-3].Fit("pol1", "0")

    canv_ratio_d_wratio.SaveAs("D_mid_fwd_ratio_vsFONLL_withRatioToFONLL.pdf")

    canv_ratio_d = ROOT.TCanvas("canv_ratio_d", "", 1400, 600)
    pad_ratio_d = [ROOT.TPad("pad_ratio_d_1", "pad_ratio_d_1", 0., 0., 0.36, 1.),
                   ROOT.TPad("pad_ratio_d_2", "pad_ratio_d_2", 0.36, 0., 0.675, 1.),
                   ROOT.TPad("pad_ratio_d_3", "pad_ratio_d_3", 0.675, 0., 1., 1.)]
    pad_ratio_d[0].SetRightMargin(0)
    pad_ratio_d[1].SetLeftMargin(0)
    pad_ratio_d[1].SetRightMargin(0)
    pad_ratio_d[2].SetLeftMargin(0)

    for irap in [2, 1, 0]:
        set_object_style(hist_ratio_stat[irap], markersize=1, color=ROOT.kAzure+4)
        set_object_style(graph_ratio_syst[irap], fillstyle=0, color=ROOT.kAzure+4)
        set_object_style(hist_ratio_d_stat[irap], markersize=1, color=ROOT.kRed+1)
        set_object_style(graph_ratio_d_syst[irap], fillstyle=0, color=ROOT.kRed+1)

    leg_ratio_b_d = ROOT.TLegend(0.04, 0.65, 0.5, 0.93)
    leg_ratio_b_d.SetBorderSize(0)
    leg_ratio_b_d.SetFillStyle(0)
    leg_ratio_b_d.SetTextSize(0.05 * size_factor[1])
    leg_ratio_b_d.AddEntry(hist_ratio_stat[-1], "#frac{B^{0} ALICE (13.6 TeV)}{B^{#plus} LHCb (13 TeV)}", "lp")
    leg_ratio_b_d.AddEntry(hist_ratio_d_stat[-1], "#frac{D^{0} ALICE (13 TeV)}{D^{0} LHCb (13 TeV)}", "lp")

    for irap in [2, 1, 0]:
        canv_ratio_d.cd()
        pad_ratio_d[irap].Draw()
        hframe = pad_ratio_d[irap].cd().DrawFrame(
            0.01, 0.3, 27., 20.,
            ";#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} mid / forward rapidity"
        )
        pad_ratio_d[irap].cd().SetLogy()
        hframe.GetYaxis().SetMoreLogLabels()
        hframe.GetYaxis().SetNoExponent()
        hframe.GetYaxis().SetTitleSize(0.05 * size_factor[irap])
        hframe.GetXaxis().SetTitleSize(0.05 * size_factor[irap])
        hframe.GetYaxis().SetLabelSize(0.05 * size_factor[irap])
        hframe.GetXaxis().SetLabelSize(0.05 * size_factor[irap])
        hframe.GetXaxis().SetTitleOffset(1.13 * 1./size_factor[irap])
        line_at_one.Draw()
        graph_ratio_syst[irap].Draw("2")
        hist_ratio_stat[irap].DrawCopy("same")
        graph_ratio_d_syst[irap].Draw("2")
        hist_ratio_d_stat[irap].DrawCopy("same")
        xmin = 0.2
        if irap > 0:
            xmin = 0.04
        else:
            lat_alice.DrawLatex(xmin+0.02, 0.88, "ALICE")
            lat_labels.DrawLatex(xmin+0.02, 0.82, "pp collisions")
        lat[irap].DrawLatex(xmin, 0.26, "ALICE |#it{y}| < 0.5")
        lat[irap].DrawLatex(xmin, 0.2, labels_lhcb[irap])
        if irap == 1:
            leg_ratio_b_d.Draw()
    canv_ratio_d.SaveAs("B_mid_fwd_ratio_vsD.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile_alice", "-a", metavar="text",
                        default="../cross_section_default_finer_pt_high_pt_w_syst.root",
                        help="input file with ALICE results", required=False)
    parser.add_argument("--infile_lhcb", "-b", metavar="text",
                        default="../lhcb/LHCb_Bplus_13TeV_rebinned.root",
                        help="input file with LHCb results", required=False)
    args = parser.parse_args()

    compute_ratio(args.infile_alice, args.infile_lhcb)

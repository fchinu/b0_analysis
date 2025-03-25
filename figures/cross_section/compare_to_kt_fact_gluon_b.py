import ROOT
import numpy as np
import sys
sys.path.append('utils')
from analysis_utils import get_n_events_from_zorro
from style_formatter import root_colors_from_matplotlib_colormap

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.12)
ROOT.gStyle.SetPadBottomMargin(0.1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

if __name__ == "__main__":
    # Get cross section
    data_file = ROOT.TFile.Open('systematics/cross_section_default_DK_pt_cuts_w_syst_fabio_fix.root')
    h_stat = data_file.Get('h_stat')
    h_stat.SetDirectory(0)
    h_syst = data_file.Get('h_syst')
    h_syst.SetDirectory(0)
    g_syst = ROOT.TGraphErrors(h_syst)
    data_file.Close()

    h_stat.Scale(1.e-6)
    h_syst.Scale(1.e-6)
    g_syst.Scale(1.e-6)

    # Get predictions
    kt_file = ROOT.TFile.Open('kt_fact/BpB0.root')
    g_pred_kt_sum = kt_file.Get('sum')
    g_pred_kt_sum.Scale(0.5)
    g_pred_kt_sum.Scale(1.e-3)
    g_pred_kt_gluon = kt_file.Get('sig2')
    g_pred_kt_gluon.Scale(0.5)
    g_pred_kt_gluon.Scale(1.e-3)
    g_pred_kt_beauty = kt_file.Get('sigm')
    g_pred_kt_beauty.Scale(0.5)
    g_pred_kt_beauty.Scale(1.e-3)
    kt_file.Close()

    pt_bins = np.append(np.asarray(g_pred_kt_sum.GetX(), 'd') - np.asarray(g_pred_kt_sum.GetEXlow(), 'd'), g_pred_kt_sum.GetX()[g_pred_kt_sum.GetN()-1] + g_pred_kt_sum.GetEXhigh()[g_pred_kt_sum.GetN()-1])

    h_pred_kt_sum = ROOT.TH1F('h_pred_kt_sum', 'h_pred_kt_sum', g_pred_kt_sum.GetN(), pt_bins)
    for i in range(1, g_pred_kt_sum.GetN()+1):
        h_pred_kt_sum.SetBinContent(i, g_pred_kt_sum.GetY()[i-1])
        h_pred_kt_sum.SetBinError(i, 1.e-10)
    
    h_pred_kt_gluon = ROOT.TH1F('h_pred_kt_gluon', 'h_pred_kt_gluon', g_pred_kt_gluon.GetN(), pt_bins)
    for i in range(1, g_pred_kt_gluon.GetN()+1):
        h_pred_kt_gluon.SetBinContent(i, g_pred_kt_gluon.GetY()[i-1])
        h_pred_kt_gluon.SetBinError(i, 1.e-10)
    
    h_pred_kt_beauty = ROOT.TH1F('h_pred_kt_beauty', 'h_pred_kt_beauty', g_pred_kt_beauty.GetN(), pt_bins)
    for i in range(1, g_pred_kt_beauty.GetN()+1):
        h_pred_kt_beauty.SetBinContent(i, g_pred_kt_beauty.GetY()[i-1])
        h_pred_kt_beauty.SetBinError(i, 1.e-10)

    # Set the style
    colors, _ = root_colors_from_matplotlib_colormap('tab10')
    h_stat.SetMarkerStyle(ROOT.kFullCircle)
    h_stat.SetMarkerSize(1.5)
    h_stat.SetMarkerColor(ROOT.kBlack)
    h_stat.SetLineColor(ROOT.kBlack)
    h_stat.SetLineWidth(2)

    g_syst.SetMarkerStyle(ROOT.kFullCircle)
    g_syst.SetMarkerSize(1.5)
    g_syst.SetMarkerColor(ROOT.kBlack)
    g_syst.SetLineColor(ROOT.kBlack)
    g_syst.SetLineWidth(2)
    g_syst.SetFillStyle(0)
    for i in range(0, g_syst.GetN()):
        g_syst.SetPointError(i, g_syst.GetErrorX(i)/2, g_syst.GetErrorY(i))

    g_pred_kt_sum.SetLineColorAlpha(colors[0], 0.5)
    g_pred_kt_sum.SetLineWidth(2)
    g_pred_kt_sum.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_kt_sum.SetFillStyle(1001)
    g_pred_kt_sum.SetFillColorAlpha(colors[0], 0.5)

    h_pred_kt_sum.SetLineColorAlpha(colors[0], 0.5)
    h_pred_kt_sum.SetLineWidth(2)
    h_pred_kt_sum.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_kt_sum.SetFillStyle(1001)
    h_pred_kt_sum.SetFillColorAlpha(colors[0], 0.5)

    g_pred_kt_gluon.SetLineColorAlpha(colors[1], 0.5)
    g_pred_kt_gluon.SetLineWidth(2)
    g_pred_kt_gluon.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_kt_gluon.SetFillStyle(1001)
    g_pred_kt_gluon.SetFillColorAlpha(colors[1], 0.5)

    h_pred_kt_gluon.SetLineColorAlpha(colors[1], 0.5)
    h_pred_kt_gluon.SetLineWidth(2)
    h_pred_kt_gluon.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_kt_gluon.SetFillStyle(1001)
    h_pred_kt_gluon.SetFillColorAlpha(colors[1], 0.5)

    g_pred_kt_beauty.SetLineColorAlpha(colors[2], 0.5)
    g_pred_kt_beauty.SetLineWidth(2)
    g_pred_kt_beauty.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_kt_beauty.SetFillStyle(1001)
    g_pred_kt_beauty.SetFillColorAlpha(colors[2], 0.5)

    h_pred_kt_beauty.SetLineColorAlpha(colors[2], 0.5)
    h_pred_kt_beauty.SetLineWidth(2)
    h_pred_kt_beauty.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_kt_beauty.SetFillStyle(1001)
    h_pred_kt_beauty.SetFillColorAlpha(colors[2], 0.5)


    # Draw
    c = ROOT.TCanvas('c', 'c', 600, 600)
    c.SetLogy()
    h_frame = c.DrawFrame(1, 2.e-3, 23.5, 3.e2, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (#mub GeV^{#minus1}#kern[0.25]{#it{c}})')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)

    g_pred_kt_sum.Draw('same E3')
    #h_pred_kt_sum.Draw('same e')
    g_pred_kt_gluon.Draw('same E3')
    #h_pred_kt_gluon.Draw('same e')
    g_pred_kt_beauty.Draw('same E3')
    #h_pred_kt_beauty.Draw('same e')
    h_stat.Draw('same p')
    g_syst.Draw('same 5')

    # Legend
    leg = ROOT.TLegend(0.4, 0.71, 0.55, 0.76)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetMargin(0.7)
    leg.AddEntry(h_stat, 'Data', 'lp')
    leg.Draw()

    # Legend
    leg_kt = ROOT.TLegend(0.65, 0.61, 0.8, 0.76)
    leg_kt.SetBorderSize(0)
    leg_kt.SetFillStyle(0)
    leg_kt.SetTextSize(0.04)
    leg_kt.SetMargin(0.7)
    leg_kt.AddEntry(g_pred_kt_gluon, 'g#rightarrowB^{0}', 'fl')
    leg_kt.AddEntry(g_pred_kt_beauty, 'b#rightarrowB^{0}', 'fl')
    leg_kt.AddEntry(g_pred_kt_sum, 'Sum', 'fl')
    leg_kt.Draw()

    # Add the text
    text_decay = ROOT.TLatex(0.41, 0.77, 'B^{0} mesons')
    text_decay.SetNDC()
    text_decay.SetTextSize(0.04)
    text_decay.SetTextFont(42)
    text_decay.Draw()

    # Add the text
    text_kt = ROOT.TLatex(0.66, 0.77, '#it{k}_{T} factorisation')
    text_kt.SetNDC()
    text_kt.SetTextSize(0.04)
    text_kt.SetTextFont(42)
    text_kt.Draw()

    text_ALICE = ROOT.TLatex(0.15, 0.88, 'ALICE Preliminary')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.055)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.15, 0.84, 'pp collisions, #sqrt{#it{s}} = 13.6 TeV')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    text_rapidity = ROOT.TLatex(0.15, 0.8, '|y| < 0.5')
    text_rapidity.SetNDC()
    text_rapidity.SetTextSize(0.04)
    text_rapidity.SetTextFont(42)
    text_rapidity.Draw()

    ROOT.gPad.RedrawAxis()

    c.SaveAs('figures/cross_section/cross_section_vs_kt_fact_gluon_b.pdf')

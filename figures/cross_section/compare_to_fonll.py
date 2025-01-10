import ROOT 
import numpy as np
import sys
sys.path.append('utils')
from analysis_utils import get_n_events_from_zorro

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.11)
ROOT.gStyle.SetPadBottomMargin(0.1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

if __name__ == '__main__':

    # Get cross section
    data_file = ROOT.TFile.Open('systematics/cross_section_default_w_syst.root')
    h_stat = data_file.Get('h_stat')
    h_stat.SetDirectory(0)
    h_syst = data_file.Get('h_syst')
    h_syst.SetDirectory(0)
    g_syst = ROOT.TGraphErrors(h_syst)
    data_file.Close()

    # Get predictions
    prediction_file = ROOT.TFile.Open('fonll/fonll_bhadron_nnpdfs_13dot6tev.root')
    g_pred_fonll = prediction_file.Get('gBhadrNNPDF30')
    prediction_file.Close()

    pt_bins = np.append(np.asarray(g_pred_fonll.GetX(), 'd') - np.asarray(g_pred_fonll.GetEXlow(), 'd'), g_pred_fonll.GetX()[g_pred_fonll.GetN()-1] + g_pred_fonll.GetEXhigh()[g_pred_fonll.GetN()-1])

    h_pred_fonll = ROOT.TH1F('h_pred', 'h_pred', g_pred_fonll.GetN(), pt_bins)
    for i in range(1, g_pred_fonll.GetN()+1):
        h_pred_fonll.SetBinContent(i, g_pred_fonll.GetY()[i-1])

    # Set the style
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

    g_pred_fonll.SetLineColorAlpha(ROOT.kRed, 0.5)
    g_pred_fonll.SetLineWidth(2)
    g_pred_fonll.SetFillStyle(1001)
    g_pred_fonll.SetFillColorAlpha(ROOT.kRed, 0.5)

    h_pred_fonll.SetLineColorAlpha(ROOT.kRed, 0.5)
    h_pred_fonll.SetLineWidth(2)
    h_pred_fonll.SetMarkerColorAlpha(ROOT.kRed, 0.5)

    # Draw
    c = ROOT.TCanvas('c', 'c', 700, 600)
    c.SetLogy()
    h_frame = c.DrawFrame(2, 5.e4, 23.5, 1.e8, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (pb #kern[-0.5]{#it{c}} / GeV)')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    g_pred_fonll.Draw('same E2')
    h_pred_fonll.Draw('same e')
    h_stat.Draw('same p')
    g_syst.Draw('same 5')

    # Legend
    leg = ROOT.TLegend(0.6, 0.65, 0.7, 0.8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetMargin(0.7)
    leg.AddEntry(h_stat, 'Data', 'lp')
    leg.AddEntry(g_pred_fonll, 'FONLL', 'fl')
    leg.Draw()

    # Add the text
    text_decay = ROOT.TLatex(0.61, 0.85, 'B^{0}#rightarrow D^{#font[122]{-}}#pi^{+}#rightarrow#pi^{#font[122]{-}}K^{+}#pi^{#font[122]{-}}#pi^{+}')
    text_decay.SetNDC()
    text_decay.SetTextSize(0.04)
    text_decay.SetTextFont(42)
    text_decay.Draw()

    text_conj = ROOT.TLatex(0.61, 0.8, 'and charge conjugate')
    text_conj.SetNDC()
    text_conj.SetTextSize(0.04)
    text_conj.SetTextFont(42)
    text_conj.Draw()

    text_ALICE = ROOT.TLatex(0.15, 0.85, 'Work in Progress')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.06)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.15, 0.8, 'pp collisions, #sqrt{#it{s}} = 13.6 TeV')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    text_rapidity = ROOT.TLatex(0.15, 0.75, '|y| < 0.5')
    text_rapidity.SetNDC()
    text_rapidity.SetTextSize(0.04)
    text_rapidity.SetTextFont(42)
    text_rapidity.Draw()
    
    ROOT.gPad.RedrawAxis()

    c.SaveAs('figures/cross_section/cross_section_vs_FONLL_23_24_full_with_bc_cuts_chebpol2_w_syst.pdf')

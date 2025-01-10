import ROOT 
import numpy as np
import sys
sys.path.append('utils')
from analysis_utils import get_n_events_from_zorro, rebin_tgraph_asymm_errors
from style_formatter import root_colors_from_matplotlib_colormap

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.11)
ROOT.gStyle.SetPadBottomMargin(0.1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

if __name__ == '__main__':

    pt_bins = [2., 4., 6., 8., 16]

    # Get cross section
    data_file = ROOT.TFile.Open('systematics/cross_section_default_w_syst.root')
    h_stat = data_file.Get('h_stat')
    h_stat.SetDirectory(0)
    h_syst = data_file.Get('h_syst')
    h_syst.SetDirectory(0)
    g_syst = ROOT.TGraphErrors(h_syst)
    data_file.Close()

    # Get CMS results
    cms_file = ROOT.TFile.Open('cms/cms_bplus_13_tev_table_2.root')
    g_stat_cms_1p45 = ROOT.TGraphAsymmErrors()
    g_stat_cms_2p1 = ROOT.TGraphAsymmErrors()
    g_syst_cms_1p45 = ROOT.TGraphAsymmErrors()
    g_syst_cms_2p1 = ROOT.TGraphAsymmErrors()

    g_stat_cms = cms_file.Get('g_stat')
    g_syst_cms = cms_file.Get('g_syst')
    for i in range(g_stat_cms.GetN()):
        if g_stat_cms.GetPointX(i) < 17:
            g_stat_cms_1p45.SetPoint(i, g_stat_cms.GetPointX(i), g_stat_cms.GetPointY(i))
            g_stat_cms_1p45.SetPointError(i, g_stat_cms.GetErrorXlow(i), g_stat_cms.GetErrorXhigh(i), g_stat_cms.GetErrorYlow(i), g_stat_cms.GetErrorYhigh(i))
            g_syst_cms_1p45.SetPoint(i, g_syst_cms.GetPointX(i), g_syst_cms.GetPointY(i))
            g_syst_cms_1p45.SetPointError(i, g_syst_cms.GetErrorXlow(i)/2, g_syst_cms.GetErrorXhigh(i)/2, g_syst_cms.GetErrorYlow(i), g_syst_cms.GetErrorYhigh(i))
        else:
            g_stat_cms_2p1.SetPoint(i, g_stat_cms.GetPointX(i), g_stat_cms.GetPointY(i))
            g_stat_cms_2p1.SetPointError(i, g_stat_cms.GetErrorXlow(i), g_stat_cms.GetErrorXhigh(i), g_stat_cms.GetErrorYlow(i), g_stat_cms.GetErrorYhigh(i))
            g_syst_cms_2p1.SetPoint(i, g_syst_cms.GetPointX(i), g_syst_cms.GetPointY(i))
            g_syst_cms_2p1.SetPointError(i, g_syst_cms.GetErrorXlow(i)/2, g_syst_cms.GetErrorXhigh(i)/2, g_syst_cms.GetErrorYlow(i), g_syst_cms.GetErrorYhigh(i))
    cms_file.Close()

    g_stat_cms_1p45.Scale(1./2.9) # divide by rapidity range
    g_stat_cms_2p1.Scale(1./4.2) # divide by rapidity range

    g_stat_cms_1p45.Scale(1.e6) # convert to pb
    g_stat_cms_2p1.Scale(1.e6) # convert to pb

    g_syst_cms_1p45.Scale(1./2.9) # divide by rapidity range
    g_syst_cms_2p1.Scale(1./4.2) # divide by rapidity range

    g_syst_cms_1p45.Scale(1.e6) # convert to pb
    g_syst_cms_2p1.Scale(1.e6) # convert to pb
    g_stat_cms_1p45.Print()

    # Set the style
    h_stat.SetMarkerStyle(ROOT.kFullCircle)
    h_stat.SetMarkerSize(1.2)
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

    g_stat_cms_1p45.SetMarkerStyle(ROOT.kOpenDiamond)
    g_stat_cms_1p45.SetMarkerSize(1.5)
    g_stat_cms_1p45.SetMarkerColor(ROOT.kRed)
    g_stat_cms_1p45.SetLineColor(ROOT.kRed)
    g_stat_cms_1p45.SetLineWidth(2)

    g_syst_cms_1p45.SetMarkerStyle(ROOT.kOpenDiamond)
    g_syst_cms_1p45.SetMarkerSize(1.5)
    g_syst_cms_1p45.SetMarkerColor(ROOT.kRed)
    g_syst_cms_1p45.SetLineColor(ROOT.kRed)
    g_syst_cms_1p45.SetLineWidth(2)
    g_syst_cms_1p45.SetFillStyle(0)

    g_stat_cms_2p1.SetMarkerStyle(ROOT.kFullDiamond)
    g_stat_cms_2p1.SetMarkerSize(1.5)
    g_stat_cms_2p1.SetMarkerColor(ROOT.kRed)
    g_stat_cms_2p1.SetLineColor(ROOT.kRed)
    g_stat_cms_2p1.SetLineWidth(2)

    g_syst_cms_2p1.SetMarkerStyle(ROOT.kFullDiamond)
    g_syst_cms_2p1.SetMarkerSize(1.5)
    g_syst_cms_2p1.SetMarkerColor(ROOT.kRed)
    g_syst_cms_2p1.SetLineColor(ROOT.kRed)
    g_syst_cms_2p1.SetLineWidth(2)
    g_syst_cms_2p1.SetFillStyle(0)

    # Draw
    c = ROOT.TCanvas('c', 'c', 700, 600)
    c.SetLogy()
    c.SetLogx()
    h_frame = c.DrawFrame(2, 5, 100, 4.e8, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (pb #it{c} / GeV)')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    g_stat_cms_1p45.Draw('same pZ')
    g_stat_cms_2p1.Draw('same pZ')
    g_syst_cms_1p45.Draw('same 5')
    g_syst_cms_2p1.Draw('same 5')
    h_stat.Draw('same p')
    g_syst.Draw('same 5')

    # Legend
    leg_alice = ROOT.TLegend(0.6, 0.68, 0.7, 0.78)
    leg_alice.SetBorderSize(0)
    leg_alice.SetFillStyle(0)
    leg_alice.SetTextSize(0.04)
    leg_alice.SetMargin(0.7)
    leg_alice.SetHeader('ALICE, #sqrt{#it{s}} = 13.6 TeV')
    leg_alice.AddEntry(h_stat, '|#it{y}| < 0.5', 'lp')
    leg_alice.Draw()

    leg_cms = ROOT.TLegend(0.14, 0.13, 0.24, 0.28)
    leg_cms.SetBorderSize(0)
    leg_cms.SetFillStyle(0)
    leg_cms.SetTextSize(0.04)
    leg_cms.SetMargin(0.7)
    leg_cms.SetHeader('CMS, #sqrt{#it{s}} = 13 TeV')
    leg_cms.AddEntry(g_stat_cms_1p45, '|#it{y}| < 1.45', 'lp')
    leg_cms.AddEntry(g_stat_cms_2p1, '|#it{y}| < 2.1', 'lp')
    leg_cms.Draw()

    # Add the text
    text_decay_alice = ROOT.TLatex(0.61, 0.85, 'B^{0}#rightarrow D^{#font[122]{-}}#pi^{+}#rightarrow #pi^{#font[122]{-}}K^{+}#pi^{#font[122]{-}}#pi^{+}')
    text_decay_alice.SetNDC()
    text_decay_alice.SetTextSize(0.04)
    text_decay_alice.SetTextFont(42)
    text_decay_alice.Draw()

    text_conj_alice = ROOT.TLatex(0.61, 0.8, 'and charge conjugate')
    text_conj_alice.SetNDC()
    text_conj_alice.SetTextSize(0.04)
    text_conj_alice.SetTextFont(42)
    text_conj_alice.Draw()

    text_decay_cms = ROOT.TLatex(0.15, 0.35, 'B^{+}#rightarrow J/#PsiK^{+}#rightarrow #mu^{#font[122]{-}}#mu^{+}K^{+}')
    text_decay_cms.SetNDC()
    text_decay_cms.SetTextSize(0.04)
    text_decay_cms.SetTextFont(42)
    text_decay_cms.Draw()

    text_conj_cms = ROOT.TLatex(0.15, 0.3, 'and charge conjugate')
    text_conj_cms.SetNDC()
    text_conj_cms.SetTextSize(0.04)
    text_conj_cms.SetTextFont(42)
    text_conj_cms.Draw()

    text_ALICE = ROOT.TLatex(0.15, 0.85, 'Work in Progress')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.06)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.15, 0.8, 'pp collisions')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    ROOT.gPad.RedrawAxis()

    c.SaveAs('figures/cross_section/cross_section_vs_CMS_23_24_w_syst.pdf')

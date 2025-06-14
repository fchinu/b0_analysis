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

if __name__ == '__main__':

    # Get cross section
    data_file = ROOT.TFile.Open('systematics/cross_section_default_DK_MC_pt_cuts_w_syst_fabio_fix_TT_vs_phi.root')
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
    tamu_file = ROOT.TFile.Open('tamu/tamu_bhadr_13dot6.root')
    g_pred_tamu = tamu_file.Get('gBhadr')
    g_pred_tamu.Scale(1.e6)
    g_pred_tamu.Scale(1.e-6)
    tamu_file.Close()

    catania_file = ROOT.TFile.Open('catania/B0_meson_13_6TeV_band.root')
    g_pred_catania = catania_file.Get('gBhadr')
    catania_file.Close()

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

    g_pred_tamu.SetLineColorAlpha(colors[3], 1)
    g_pred_tamu.SetLineWidth(3)
    g_pred_tamu.SetFillStyle(1001)
    g_pred_tamu.SetLineStyle(9)
    g_pred_tamu.SetFillColorAlpha(colors[3], 0.5)

    g_pred_catania.SetLineColorAlpha(colors[5], 0.5)
    g_pred_catania.SetLineWidth(0)
    g_pred_catania.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_catania.SetFillStyle(1001)
    g_pred_catania.SetFillColorAlpha(colors[5], 0.7)

    # Draw
    c = ROOT.TCanvas('c', 'c', 900, 600)
    pad_cross_sec = ROOT.TPad('pad_cross_sec', 'pad_cross_sec', 0, 0, 0.5, 1)
    pad_cross_sec.Draw()
    pad_cross_sec.cd()
    pad_cross_sec.SetLogy()
    h_frame = pad_cross_sec.DrawFrame(1, 2.e-2, 23.5, 3.e2, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (#mub GeV^{#minus1}#kern[0.25]{#it{c}})')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    g_pred_tamu.Draw('same L')
    g_pred_catania.Draw('same E2')
    h_stat.Draw('same p')
    g_syst.Draw('same 5')

    # Legend
    leg = ROOT.TLegend(0.4, 0.635, 0.55, 0.755)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetMargin(0.7)
    leg.AddEntry(h_stat, 'Data', 'lp')
    leg.AddEntry(g_pred_catania, 'Catania', 'f')
    leg.AddEntry(g_pred_tamu, 'TAMU', 'l')
    leg.Draw()

    # Add the text
    text_decay = ROOT.TLatex(0.41, 0.765, 'B^{0} mesons')
    text_decay.SetNDC()
    text_decay.SetTextSize(0.04)
    text_decay.SetTextFont(42)
    text_decay.Draw()

    text_ALICE = ROOT.TLatex(0.15, 0.885, 'ALICE')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.06)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.15, 0.855, 'pp collisions, #sqrt{#it{s}} = 13.6 TeV')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    text_rapidity = ROOT.TLatex(0.15, 0.82, '|y| < 0.5')
    text_rapidity.SetNDC()
    text_rapidity.SetTextSize(0.04)
    text_rapidity.SetTextFont(42)
    text_rapidity.Draw()

    text_lumi = ROOT.TLatex(0.15, 0.785, '#font[132]{#it{L}}_{int} = 43 pb^{#minus1}')
    text_lumi.SetNDC()
    text_lumi.SetTextSize(0.04) 
    text_lumi.SetTextFont(42)
    text_lumi.Draw()

    ROOT.gPad.RedrawAxis()


    c.cd()
    pad_ratios_axis_label = ROOT.TPad('pad_ratios_axis_label', 'pad_ratios_axis_label', 0.5, 0, 1, 1)
    pad_ratios_axis_label.SetLeftMargin(0.15)
    pad_ratios_axis_label.Draw()
    pad_ratios_axis_label.cd()
    h_frame_axis_label = pad_ratios_axis_label.DrawFrame(1, 0.123, 23.5, 0.9876, ';#it{p}_{T} (GeV/#it{c});')
    h_frame_axis_label.GetXaxis().SetTitleOffset(1.1)
    h_frame_axis_label.GetYaxis().SetTitleOffset(1.3)
    h_frame_axis_label.GetXaxis().SetTitleSize(0.04)
    h_frame_axis_label.GetYaxis().SetTitleSize(0.04)
    h_frame_axis_label.GetXaxis().SetLabelSize(0.04)
    h_frame_axis_label.GetYaxis().SetLabelSize(0.04)

    c.cd()

    pad_ratios = ROOT.TPad('pad_ratios', 'pad_ratios', 0.5, 0.1, 1, 0.95)
    pad_ratios.Draw()
    pad_ratios.cd()

    pad_ratio_catania = ROOT.TPad('pad_ratio_catania', 'pad_ratio_catania', 0, 0.5, 1, 1)
    pad_ratio_catania.Draw()
    pad_ratio_catania.cd()
    pad_ratio_catania.SetBottomMargin(0)
    pad_ratio_catania.SetLeftMargin(0.15)
    pad_ratio_catania.SetTopMargin(0)
    h_frame_ratio_catania = pad_ratio_catania.DrawFrame(1, 0.2, 23.5,2.995, ';#it{p}_{T} (GeV/#it{c});#frac{Data}{Catania}')
    h_frame_ratio_catania.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_catania.GetYaxis().SetTitleOffset(0.9)
    h_frame_ratio_catania.GetYaxis().SetDecimals(True)
    h_frame_ratio_catania.GetYaxis().CenterTitle(True)
    h_frame_ratio_catania.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_catania.GetYaxis().SetTitleSize(0.08)
    h_frame_ratio_catania.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_catania.GetYaxis().SetLabelSize(0.07)
    h_frame_ratio_catania.GetYaxis().SetNdivisions(507)

    h_catania_int = h_stat.Clone()
    h_catania_int.Reset()
    catania_step = 0.2
    for i in range(0, g_pred_catania.GetN()):
        h_catania_int.Fill(g_pred_catania.GetX()[i], g_pred_catania.GetY()[i] * catania_step / h_catania_int.GetBinWidth(h_catania_int.FindBin(g_pred_catania.GetX()[i])))
    h_catania_int.SetBinContent(h_catania_int.GetNbinsX()-1, 1.e-12)

    for i in range(1, h_catania_int.GetNbinsX()+1):
        err = 0
        for j in range(0, g_pred_catania.GetN()):
            if g_pred_catania.GetX()[j] >= h_catania_int.GetXaxis().GetBinLowEdge(i) and g_pred_catania.GetX()[j] < h_catania_int.GetXaxis().GetBinUpEdge(i):
                err += g_pred_catania.GetErrorYhigh(j)

        h_catania_int.SetBinError(i, catania_step * err / h_catania_int.GetBinWidth(i))
    h_catania_int.SetBinError(h_catania_int.GetNbinsX()-1, 0)

    # pad_cross_sec.cd()
    # h_catania_int.Draw('same e2')
    # pad_ratio_catania.cd()

    pt_bins = h_catania_int.GetXaxis().GetXbins()

    g_catania_unc = ROOT.TGraphAsymmErrors(h_catania_int.Clone('g_catania_unc'))
    for i in range(0, g_catania_unc.GetN()-1):
        g_catania_unc.SetPoint(i, g_catania_unc.GetX()[i], 1)
        g_catania_unc.SetPointError(i, g_catania_unc.GetErrorXlow(i), g_catania_unc.GetErrorXhigh(i), g_catania_unc.GetErrorYlow(i)/h_catania_int.GetBinContent(i+1), g_catania_unc.GetErrorYhigh(i)/h_catania_int.GetBinContent(i+1))

    g_catania_unc.SetMarkerSize(0)
    g_catania_unc.SetLineColorAlpha(colors[5], 1)
    g_catania_unc.SetLineWidth(0)
    g_catania_unc.SetFillStyle(1001)
    # g_catania_unc.SetLineStyle(9)
    g_catania_unc.SetFillColorAlpha(colors[5], 0.5)

    g_catania_unc.Draw("same 5")

    line_one_catania = ROOT.TLine(1, 1, 23.5, 1)
    line_one_catania.SetLineStyle(2)
    line_one_catania.SetLineColor(ROOT.kBlack)
    line_one_catania.Draw("same")

    h_ratio_data_catania_stat = h_stat.Clone('h_ratio_data_catania_stat')
    for i in range(1, h_ratio_data_catania_stat.GetNbinsX()):
        h_ratio_data_catania_stat.SetBinContent(i, h_stat.GetBinContent(i)/h_catania_int.GetBinContent(i))
        h_ratio_data_catania_stat.SetBinError(i, h_stat.GetBinError(i)/h_catania_int.GetBinContent(i))
        h_ratio_data_catania_stat.Draw('same p')

    g_ratio_data_catania_syst = g_syst.Clone('g_ratio_data_catania_syst')
    for i in range(0, g_ratio_data_catania_syst.GetN()-1):
        g_ratio_data_catania_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/h_catania_int.GetBinContent(i+1))
        g_ratio_data_catania_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/h_catania_int.GetBinContent(i+1))

    g_ratio_data_catania_syst.Draw('same 5') 

    pad_ratio_catania.RedrawAxis()

    pad_ratios.cd()

    pad_ratio_tamu = ROOT.TPad('pad_ratio_tamu', 'pad_ratio_tamu', 0, 0., 1, 0.5)
    pad_ratio_tamu.Draw()
    pad_ratio_tamu.cd()
    pad_ratio_tamu.SetBottomMargin(0)
    pad_ratio_tamu.SetLeftMargin(0.15)
    pad_ratio_tamu.SetTopMargin(0)
    h_frame_ratio_tamu = pad_ratio_tamu.DrawFrame(1, 0.2, 23.5,2.995, ';#it{p}_{T} (GeV/#it{c});#frac{Data}{TAMU}')
    h_frame_ratio_tamu.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_tamu.GetYaxis().SetTitleOffset(0.8)
    h_frame_ratio_tamu.GetYaxis().SetDecimals(True)
    h_frame_ratio_tamu.GetYaxis().CenterTitle(True)
    h_frame_ratio_tamu.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_tamu.GetYaxis().SetTitleSize(0.08)
    h_frame_ratio_tamu.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_tamu.GetYaxis().SetLabelSize(0.08)
    h_frame_ratio_tamu.GetYaxis().SetNdivisions(507)

    h_tamu_int = h_stat.Clone()
    h_tamu_int.Reset()
    tamu_step = 0.2
    for i in range(0, g_pred_tamu.GetN()):
        h_tamu_int.Fill(g_pred_tamu.GetX()[i], g_pred_tamu.GetY()[i] * tamu_step / h_tamu_int.GetBinWidth(h_tamu_int.FindBin(g_pred_tamu.GetX()[i])))
    h_tamu_int.SetBinContent(h_tamu_int.GetNbinsX(), 1.e-12)

    for i in range(1, h_tamu_int.GetNbinsX()+1):
        h_tamu_int.SetBinError(i, 0)

    pt_bins = h_tamu_int.GetXaxis().GetXbins()
    h_tamu_unc = ROOT.TH1D('h_tamu_unc', 'h_tamu_unc', g_pred_tamu.GetN(), 0, 20.)
    for i in range(1, h_tamu_unc.GetNbinsX()+1):
        h_tamu_unc.SetBinContent(i, 1)
        h_tamu_unc.SetBinError(i, 0.)

    h_tamu_unc.SetMarkerSize(0)
    h_tamu_unc.SetLineColorAlpha(colors[3], 1)
    h_tamu_unc.SetLineWidth(3)
    h_tamu_unc.SetFillStyle(0)
    h_tamu_unc.SetLineStyle(9)
    h_tamu_unc.SetFillColorAlpha(colors[3], 0.5)

    line_one_tamu = ROOT.TLine(1, 1, 23.5, 1)
    line_one_tamu.SetLineStyle(2)
    line_one_tamu.SetLineColor(ROOT.kBlack)
    line_one_tamu.Draw("same")

    h_tamu_unc.Draw("same l")

    h_ratio_data_tamu_stat = h_stat.Clone('h_ratio_data_tamu_stat')
    h_ratio_data_tamu_stat.Divide(h_tamu_int)
    h_ratio_data_tamu_stat.Draw('same p')

    g_ratio_data_tamu_syst = g_syst.Clone('g_ratio_data_tamu_syst')
    for i in range(0, g_ratio_data_tamu_syst.GetN()):
        g_ratio_data_tamu_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/h_tamu_int.GetBinContent(i+1))
        g_ratio_data_tamu_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/h_tamu_int.GetBinContent(i+1))

    g_ratio_data_tamu_syst.Draw('same 5') 

    pad_ratio_tamu.RedrawAxis()

    c.SaveAs('figures/cross_section/cross_section_vs_pheno_TT_dk_mc.pdf')

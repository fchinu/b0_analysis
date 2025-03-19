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
    fonll_file = ROOT.TFile.Open('fonll/fonll_bhadron_nnpdfs_13dot6tev.root')
    g_pred_fonll = fonll_file.Get('gBhadrNNPDF30')
    fonll_file.Close()
    g_pred_fonll.Scale(1.e-6)

    gmvfns_file = ROOT.TFile.Open('gmvfns/gmvfns_pred_13dot6.root')
    g_pred_gmvfns = gmvfns_file.Get('gBhadr')
    g_pred_gmvfns.Scale(1.e3)
    g_pred_gmvfns.Scale(1.e-6)
    gmvfns_file.Close()

    gmvfns_mt_sacot_file = ROOT.TFile.Open('sacot_mt/sacot_mt_pred_13dot6.root')
    g_pred_gmvfns_mt_sacot = gmvfns_mt_sacot_file.Get('gBhadr')
    g_pred_gmvfns_mt_sacot.Scale(1.e9)
    g_pred_gmvfns_mt_sacot.Scale(1.e-6)
    gmvfns_mt_sacot_file.Close()

    kt_file = ROOT.TFile.Open('kt_fact/BpB0_binning.root')
    g_pred_kt = kt_file.Get('sum')
    g_pred_kt.Scale(0.5)
    g_pred_kt.Scale(1.e-3)
    kt_file.Close()

    tamu_file = ROOT.TFile.Open('tamu/tamu_bhadr_13dot6.root')
    g_pred_tamu = tamu_file.Get('gBhadr')
    g_pred_tamu.Scale(1.e6)
    g_pred_tamu.Scale(1.e-6)
    tamu_file.Close()

    pt_bins = np.append(np.asarray(g_pred_fonll.GetX(), 'd') - np.asarray(g_pred_fonll.GetEXlow(), 'd'), g_pred_fonll.GetX()[g_pred_fonll.GetN()-1] + g_pred_fonll.GetEXhigh()[g_pred_fonll.GetN()-1])

    h_pred_fonll = ROOT.TH1F('h_pred_fonll', 'h_pred_fonll', g_pred_fonll.GetN(), pt_bins)
    for i in range(1, g_pred_fonll.GetN()+1):
        h_pred_fonll.SetBinContent(i, g_pred_fonll.GetY()[i-1])
        h_pred_fonll.SetBinError(i, 1.e-10)

    h_pred_gmvfns = ROOT.TH1F('h_pred_gmvfns', 'h_pred_gmvfns', g_pred_gmvfns.GetN(), pt_bins)
    for i in range(1, g_pred_gmvfns.GetN()+1):
        h_pred_gmvfns.SetBinContent(i, g_pred_gmvfns.GetY()[i-1])
        h_pred_gmvfns.SetBinError(i, 1.e-10)

    h_pred_gmvfns_mt_sacot = ROOT.TH1F('h_pred_gmvfns_mt_sacot', 'h_pred_gmvfns_mt_sacot', g_pred_gmvfns_mt_sacot.GetN(), pt_bins)
    for i in range(1, g_pred_gmvfns_mt_sacot.GetN()+1):
        h_pred_gmvfns_mt_sacot.SetBinContent(i, g_pred_gmvfns_mt_sacot.GetY()[i-1])
        h_pred_gmvfns_mt_sacot.SetBinError(i, 1.e-10)

    h_pred_kt = ROOT.TH1F('h_pred_kt', 'h_pred_kt', g_pred_kt.GetN(), pt_bins)
    for i in range(1, g_pred_kt.GetN()+1):
        h_pred_kt.SetBinContent(i, g_pred_kt.GetY()[i-1])
        h_pred_kt.SetBinError(i, 1.e-10)

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

    g_pred_fonll.SetLineColorAlpha(colors[0], 0.5)
    g_pred_fonll.SetLineWidth(2)
    g_pred_fonll.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_fonll.SetFillStyle(1001)
    g_pred_fonll.SetFillColorAlpha(colors[0], 0.7)

    h_pred_fonll.SetLineColorAlpha(colors[0], 0.5)
    h_pred_fonll.SetLineWidth(2)
    h_pred_fonll.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_fonll.SetFillColorAlpha(colors[0], 0.7)
    
    g_pred_gmvfns.SetLineColorAlpha(colors[1], 0.5)
    g_pred_gmvfns.SetLineWidth(2)
    g_pred_gmvfns.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_gmvfns.SetFillStyle(1001)
    g_pred_gmvfns.SetFillColorAlpha(colors[1], 0.7)

    h_pred_gmvfns.SetLineColorAlpha(colors[1], 0.5)
    h_pred_gmvfns.SetLineWidth(2)
    h_pred_gmvfns.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_gmvfns.SetFillColorAlpha(colors[1], 0.7)

    g_pred_gmvfns_mt_sacot.SetLineColorAlpha(colors[2], 0.5)
    g_pred_gmvfns_mt_sacot.SetLineWidth(2)
    g_pred_gmvfns_mt_sacot.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_gmvfns_mt_sacot.SetFillStyle(1001)
    g_pred_gmvfns_mt_sacot.SetFillColorAlpha(colors[2], 0.35)

    h_pred_gmvfns_mt_sacot.SetLineColorAlpha(colors[2], 0.5)
    h_pred_gmvfns_mt_sacot.SetLineWidth(2)
    h_pred_gmvfns_mt_sacot.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_gmvfns_mt_sacot.SetFillColorAlpha(colors[2], 0.35)

    g_pred_kt.SetLineColorAlpha(colors[4], 0.5)
    g_pred_kt.SetLineWidth(2)
    g_pred_kt.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    g_pred_kt.SetFillStyle(1001)
    g_pred_kt.SetFillColorAlpha(colors[4], 0.35)

    h_pred_kt.SetLineColorAlpha(colors[4], 0.5)
    h_pred_kt.SetLineWidth(2)
    h_pred_kt.SetMarkerColorAlpha(ROOT.kBlack, 0.)
    h_pred_kt.SetFillColorAlpha(colors[4], 0.35)

    g_pred_tamu.SetLineColorAlpha(colors[3], 1)
    g_pred_tamu.SetLineWidth(3)
    g_pred_tamu.SetFillStyle(1001)
    g_pred_tamu.SetLineStyle(9)
    g_pred_tamu.SetFillColorAlpha(colors[3], 0.5)

    # Draw
    c = ROOT.TCanvas('c', 'c', 900, 600)
    pad_cross_sec = ROOT.TPad('pad_cross_sec', 'pad_cross_sec', 0, 0, 0.5, 1)
    pad_cross_sec.Draw()
    pad_cross_sec.cd()
    pad_cross_sec.SetLogy()
    h_frame = pad_cross_sec.DrawFrame(1, 2.e-2, 23.5, 3.e2, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (#mub #kern[-0.5]{#it{c}} / GeV)')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    g_pred_kt.Draw('same E2')
    h_pred_kt.Draw('same e')
    g_pred_gmvfns_mt_sacot.Draw('same E2')
    h_pred_gmvfns_mt_sacot.Draw('same e')
    g_pred_fonll.Draw('same E2')
    h_pred_fonll.Draw('same e')
    g_pred_gmvfns.Draw('same E2')
    h_pred_gmvfns.Draw('same e')
    g_pred_tamu.Draw('same L')
    h_stat.Draw('same p')
    g_syst.Draw('same 5')

    # Legend
    leg = ROOT.TLegend(0.4, 0.56, 0.55, 0.77)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetMargin(0.7)
    leg.AddEntry(h_stat, 'Data', 'lp')
    leg.AddEntry(g_pred_fonll, 'FONLL', 'fl')
    leg.AddEntry(g_pred_gmvfns, 'GM-VFNS (mod-#mu_{#lower[-0.2]{R,F}})', 'fl')
    leg.AddEntry(g_pred_gmvfns_mt_sacot, 'GM-VFNS, SACOT#kern[0.3]{#it{m}_{T}}', 'fl')
    leg.AddEntry(g_pred_kt, '#it{k}_{T} fact.', 'fl')
    leg.AddEntry(g_pred_tamu, 'TAMU', 'l')
    leg.Draw()

    # Add the text
    text_decay = ROOT.TLatex(0.41, 0.78, 'B^{0} mesons')
    text_decay.SetNDC()
    text_decay.SetTextSize(0.04)
    text_decay.SetTextFont(42)
    text_decay.Draw()

    text_ALICE = ROOT.TLatex(0.15, 0.88, 'ALICE Preliminary')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.06)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.15, 0.845, 'pp collisions, #sqrt{#it{s}} = 13.6 TeV')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    text_rapidity = ROOT.TLatex(0.15, 0.81, '|y| < 0.5')
    text_rapidity.SetNDC()
    text_rapidity.SetTextSize(0.04)
    text_rapidity.SetTextFont(42)
    text_rapidity.Draw()

    ROOT.gPad.RedrawAxis()


    c.cd()
    pad_ratios_axis_label = ROOT.TPad('pad_ratios_axis_label', 'pad_ratios_axis_label', 0.5, 0, 1, 1)
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

    pad_ratio_fonll = ROOT.TPad('pad_ratio_fonll', 'pad_ratio_fonll', 0, 0.8, 1, 1)
    pad_ratio_fonll.Draw()
    pad_ratio_fonll.cd()
    pad_ratio_fonll.SetBottomMargin(0)
    pad_ratio_fonll.SetTopMargin(0)
    h_frame_ratio_fonll = pad_ratio_fonll.DrawFrame(1, 0.2, 23.5,2.495, ';#it{p}_{T} (GeV/#it{c});#frac{Data}{FONLL}')
    h_frame_ratio_fonll.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_fonll.GetYaxis().SetTitleOffset(0.35)
    h_frame_ratio_fonll.GetYaxis().CenterTitle(True)
    h_frame_ratio_fonll.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_fonll.GetYaxis().SetTitleSize(0.16)
    h_frame_ratio_fonll.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_fonll.GetYaxis().SetLabelSize(0.16)
    h_frame_ratio_fonll.GetYaxis().SetNdivisions(505)

    g_fonll_unc = g_pred_fonll.Clone('g_fonll_unc')
    for i in range(0, g_fonll_unc.GetN()):
        g_fonll_unc.SetPoint(i, g_pred_fonll.GetX()[i], 1)
        g_fonll_unc.SetPointError(i, g_pred_fonll.GetErrorXlow(i), g_pred_fonll.GetErrorXhigh(i), g_pred_fonll.GetErrorYlow(i)/g_pred_fonll.GetY()[i], g_pred_fonll.GetErrorYhigh(i)/g_pred_fonll.GetY()[i])

    g_fonll_unc.Draw("same 2")

    line_one_fonll = ROOT.TLine(1, 1, 23.5, 1)
    line_one_fonll.SetLineStyle(2)
    line_one_fonll.SetLineColor(ROOT.kBlack)
    line_one_fonll.Draw("same")

    h_ratio_data_fonll_stat = h_stat.Clone('h_ratio_data_fonll_stat')
    h_ratio_data_fonll_stat.Divide(h_pred_fonll)
    h_ratio_data_fonll_stat.Draw('same p')

    g_ratio_data_fonll_syst = g_syst.Clone('g_ratio_data_fonll_syst')
    for i in range(0, g_ratio_data_fonll_syst.GetN()):
        g_ratio_data_fonll_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/g_pred_fonll.GetY()[i])
        g_ratio_data_fonll_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/g_pred_fonll.GetY()[i])

    g_ratio_data_fonll_syst.Draw('same 5') 

    pad_ratio_fonll.RedrawAxis()

    pad_ratios.cd()

    pad_ratio_gmvfns = ROOT.TPad('pad_ratio_gmvfns', 'pad_ratio_gmvfns', 0, 0.6, 1, 0.8)
    pad_ratio_gmvfns.Draw()
    pad_ratio_gmvfns.cd()
    pad_ratio_gmvfns.SetBottomMargin(0)
    pad_ratio_gmvfns.SetTopMargin(0)
    h_frame_ratio_gmvfns = pad_ratio_gmvfns.DrawFrame(1, 0.2, 23.5,2.495, ';#it{p}_{T} (GeV/#it{c});#lower[0.1]{#frac{Data}{mod-#mu_{#lower[-0.2]{R,F}}}}')
    h_frame_ratio_gmvfns.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_gmvfns.GetYaxis().SetTitleOffset(0.35)
    h_frame_ratio_gmvfns.GetYaxis().CenterTitle(True)
    h_frame_ratio_gmvfns.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_gmvfns.GetYaxis().SetTitleSize(0.16)
    h_frame_ratio_gmvfns.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_gmvfns.GetYaxis().SetLabelSize(0.16)
    h_frame_ratio_gmvfns.GetYaxis().SetNdivisions(505)

    g_gmvfns_unc = g_pred_gmvfns.Clone('g_gmvfns_unc')
    for i in range(0, g_gmvfns_unc.GetN()):
        g_gmvfns_unc.SetPoint(i, g_pred_gmvfns.GetX()[i], 1)
        g_gmvfns_unc.SetPointError(i, g_pred_gmvfns.GetErrorXlow(i), g_pred_gmvfns.GetErrorXhigh(i), g_pred_gmvfns.GetErrorYlow(i)/g_pred_gmvfns.GetY()[i], g_pred_gmvfns.GetErrorYhigh(i)/g_pred_gmvfns.GetY()[i])

    g_gmvfns_unc.Draw("same 2")

    line_one_gmvfns = ROOT.TLine(1, 1, 23.5, 1)
    line_one_gmvfns.SetLineStyle(2)
    line_one_gmvfns.SetLineColor(ROOT.kBlack)
    line_one_gmvfns.Draw("same")

    h_ratio_data_gmvfns_stat = h_stat.Clone('h_ratio_data_gmvfns_stat')
    h_ratio_data_gmvfns_stat.Divide(h_pred_gmvfns)
    h_ratio_data_gmvfns_stat.Draw('same p')

    g_ratio_data_gmvfns_syst = g_syst.Clone('g_ratio_data_gmvfns_syst')
    for i in range(0, g_ratio_data_gmvfns_syst.GetN()):
        g_ratio_data_gmvfns_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/g_pred_gmvfns.GetY()[i])
        g_ratio_data_gmvfns_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/g_pred_gmvfns.GetY()[i])

    g_ratio_data_gmvfns_syst.Draw('same 5')

    pad_ratio_gmvfns.RedrawAxis()

    pad_ratios.cd()

    pad_ratio_gmvfns_mt_sacot = ROOT.TPad('pad_ratio_gmvfns_mt_sacot', 'pad_ratio_gmvfns_mt_sacot', 0, 0.4, 1, 0.6)
    pad_ratio_gmvfns_mt_sacot.Draw()
    pad_ratio_gmvfns_mt_sacot.cd()
    pad_ratio_gmvfns_mt_sacot.SetBottomMargin(0)
    pad_ratio_gmvfns_mt_sacot.SetTopMargin(0)
    h_frame_ratio_gmvfns_mt_sacot = pad_ratio_gmvfns_mt_sacot.DrawFrame(1, 0.2, 23.5,2.495, ';#it{p}_{T} (GeV/#it{c});#lower[0.1]{#frac{Data}{SACOT#kern[0.3]{#it{m}_{T}}}}')
    h_frame_ratio_gmvfns_mt_sacot.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_gmvfns_mt_sacot.GetYaxis().SetTitleOffset(0.35)
    h_frame_ratio_gmvfns_mt_sacot.GetYaxis().CenterTitle(True)
    h_frame_ratio_gmvfns_mt_sacot.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_gmvfns_mt_sacot.GetYaxis().SetTitleSize(0.16)
    h_frame_ratio_gmvfns_mt_sacot.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_gmvfns_mt_sacot.GetYaxis().SetLabelSize(0.16)
    h_frame_ratio_gmvfns_mt_sacot.GetYaxis().SetNdivisions(505)

    g_gmvfns_mt_sacot_unc = g_pred_gmvfns_mt_sacot.Clone('g_gmvfns_mt_sacot_unc')
    for i in range(0, g_gmvfns_mt_sacot_unc.GetN()):
        g_gmvfns_mt_sacot_unc.SetPoint(i, g_pred_gmvfns_mt_sacot.GetX()[i], 1)
        g_gmvfns_mt_sacot_unc.SetPointError(i, g_pred_gmvfns_mt_sacot.GetErrorXlow(i), g_pred_gmvfns_mt_sacot.GetErrorXhigh(i), g_pred_gmvfns_mt_sacot.GetErrorYlow(i)/g_pred_gmvfns_mt_sacot.GetY()[i], g_pred_gmvfns_mt_sacot.GetErrorYhigh(i)/g_pred_gmvfns_mt_sacot.GetY()[i])

    g_gmvfns_mt_sacot_unc.Draw("same 2")

    line_one_gmvfns_mt_sacot = ROOT.TLine(1, 1, 23.5, 1)
    line_one_gmvfns_mt_sacot.SetLineStyle(2)
    line_one_gmvfns_mt_sacot.SetLineColor(ROOT.kBlack)
    line_one_gmvfns_mt_sacot.Draw("same")

    h_ratio_data_gmvfns_mt_sacot_stat = h_stat.Clone('h_ratio_data_gmvfns_mt_sacot_stat')
    h_ratio_data_gmvfns_mt_sacot_stat.Divide(h_pred_gmvfns_mt_sacot)
    h_ratio_data_gmvfns_mt_sacot_stat.Draw('same p')

    g_ratio_data_gmvfns_mt_sacot_syst = g_syst.Clone('g_ratio_data_gmvfns_mt_sacot_syst')
    for i in range(0, g_ratio_data_gmvfns_mt_sacot_syst.GetN()):
        g_ratio_data_gmvfns_mt_sacot_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/g_pred_gmvfns_mt_sacot.GetY()[i])
        g_ratio_data_gmvfns_mt_sacot_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/g_pred_gmvfns_mt_sacot.GetY()[i])

    g_ratio_data_gmvfns_mt_sacot_syst.Draw('same 5') 

    pad_ratio_gmvfns_mt_sacot.RedrawAxis()

    pad_ratios.cd()

    pad_ratio_kt = ROOT.TPad('pad_ratio_kt', 'pad_ratio_kt', 0, 0.2, 1, 0.4)
    pad_ratio_kt.Draw()
    pad_ratio_kt.cd()
    pad_ratio_kt.SetBottomMargin(0)
    pad_ratio_kt.SetTopMargin(0)
    h_frame_ratio_kt = pad_ratio_kt.DrawFrame(1, 0.2, 23.5,2.495, ';#it{p}_{T} (GeV/#it{c});#lower[0.09]{#frac{Data}{#it{k}_{T} fact.}}')
    h_frame_ratio_kt.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_kt.GetYaxis().SetTitleOffset(0.35)
    h_frame_ratio_kt.GetYaxis().CenterTitle(True)
    h_frame_ratio_kt.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_kt.GetYaxis().SetTitleSize(0.16)
    h_frame_ratio_kt.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_kt.GetYaxis().SetLabelSize(0.16)
    h_frame_ratio_kt.GetYaxis().SetNdivisions(505)

    g_kt_unc = g_pred_kt.Clone('g_kt_unc')
    for i in range(0, g_kt_unc.GetN()):
        g_kt_unc.SetPoint(i, g_pred_kt.GetX()[i], 1)
        g_kt_unc.SetPointError(i, g_pred_kt.GetErrorXlow(i), g_pred_kt.GetErrorXhigh(i), g_pred_kt.GetErrorYlow(i)/g_pred_kt.GetY()[i], g_pred_kt.GetErrorYhigh(i)/g_pred_kt.GetY()[i])

    g_kt_unc.Draw("same 2")

    line_one_kt = ROOT.TLine(1, 1, 23.5, 1)
    line_one_kt.SetLineStyle(2)
    line_one_kt.SetLineColor(ROOT.kBlack)
    line_one_kt.Draw("same")

    h_ratio_data_kt_stat = h_stat.Clone('h_ratio_data_kt_stat')
    h_ratio_data_kt_stat.Divide(h_pred_kt)
    h_ratio_data_kt_stat.Draw('same p')

    g_ratio_data_kt_syst = g_syst.Clone('g_ratio_data_kt_syst')
    for i in range(0, g_ratio_data_kt_syst.GetN()):
        g_ratio_data_kt_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/g_pred_kt.GetY()[i])
        g_ratio_data_kt_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/g_pred_kt.GetY()[i])

    g_ratio_data_kt_syst.Draw('same 5') 

    pad_ratio_kt.RedrawAxis()

    pad_ratios.cd()

    pad_ratio_tamu = ROOT.TPad('pad_ratio_tamu', 'pad_ratio_tamu', 0, 0., 1, 0.2)
    pad_ratio_tamu.Draw()
    pad_ratio_tamu.cd()
    pad_ratio_tamu.SetBottomMargin(0)
    pad_ratio_tamu.SetTopMargin(0)
    h_frame_ratio_tamu = pad_ratio_tamu.DrawFrame(1, 0.2, 23.5,2.495, ';#it{p}_{T} (GeV/#it{c});#frac{Data}{TAMU}')
    h_frame_ratio_tamu.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_tamu.GetYaxis().SetTitleOffset(0.35)
    h_frame_ratio_tamu.GetYaxis().CenterTitle(True)
    h_frame_ratio_tamu.GetXaxis().SetTitleSize(0.12)
    h_frame_ratio_tamu.GetYaxis().SetTitleSize(0.16)
    h_frame_ratio_tamu.GetXaxis().SetLabelSize(0.12)
    h_frame_ratio_tamu.GetYaxis().SetLabelSize(0.16)
    h_frame_ratio_tamu.GetYaxis().SetNdivisions(505)

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

    c.SaveAs('figures/cross_section/cross_section_vs_pred.pdf')

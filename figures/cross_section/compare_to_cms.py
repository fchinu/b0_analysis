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
    data_file = ROOT.TFile.Open('systematics/cross_section_default_DK_pt_cuts_w_syst_fabio_fix.root')
    h_stat = data_file.Get('h_stat')
    h_stat.SetDirectory(0)
    h_syst = data_file.Get('h_syst')
    h_syst.SetDirectory(0)
    g_syst = ROOT.TGraphErrors(h_syst)
    data_file.Close()

    h_stat.Scale(1.e-6) # convert to µb
    g_syst.Scale(1.e-6) # convert to µb

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
            g_stat_cms_2p1.SetPoint(i-2, g_stat_cms.GetPointX(i), g_stat_cms.GetPointY(i))
            g_stat_cms_2p1.SetPointError(i-2, g_stat_cms.GetErrorXlow(i), g_stat_cms.GetErrorXhigh(i), g_stat_cms.GetErrorYlow(i), g_stat_cms.GetErrorYhigh(i))
            g_syst_cms_2p1.SetPoint(i-2, g_syst_cms.GetPointX(i), g_syst_cms.GetPointY(i))
            g_syst_cms_2p1.SetPointError(i-2, g_syst_cms.GetErrorXlow(i)/2, g_syst_cms.GetErrorXhigh(i)/2, g_syst_cms.GetErrorYlow(i), g_syst_cms.GetErrorYhigh(i))
    cms_file.Close()

    g_stat_cms_1p45.Scale(1./2.9) # divide by rapidity range
    g_stat_cms_2p1.Scale(1./4.2) # divide by rapidity range

    g_syst_cms_1p45.Scale(1./2.9) # divide by rapidity range
    g_syst_cms_2p1.Scale(1./4.2) # divide by rapidity range

    g_stat_cms_1p45.Print()


    # Get predictions
    fonll_file = ROOT.TFile.Open('fonll/fonll_bhadron_nnpdfs_13dot6tev.root')
    g_pred_fonll_alice = fonll_file.Get('gBhadrNNPDF30')
    fonll_file.Close()

    g_pred_fonll_alice.Scale(1.e-6) # convert to µb

    pt_bins_fonll = np.append(np.asarray(g_pred_fonll_alice.GetX(), 'd') - np.asarray(g_pred_fonll_alice.GetEXlow(), 'd'), g_pred_fonll_alice.GetX()[g_pred_fonll_alice.GetN()-1] + g_pred_fonll_alice.GetEXhigh()[g_pred_fonll_alice.GetN()-1])

    h_pred_fonll_alice = ROOT.TH1F('h_pred_fonll', 'h_pred_fonll', g_pred_fonll_alice.GetN(), pt_bins_fonll)
    for i in range(1, g_pred_fonll_alice.GetN()+1):
        h_pred_fonll_alice.SetBinContent(i, g_pred_fonll_alice.GetY()[i-1])
        h_pred_fonll_alice.SetBinError(i, 1.e-12)

    # Get predictions
    fonll_file = ROOT.TFile.Open('fonll/fonll_bhadron_nnpdfs_13tev_y1p45.root')
    g_pred_fonll_cms_1p45 = fonll_file.Get('gBhadrNNPDF30')
    fonll_file.Close()

    g_pred_fonll_cms_1p45.Scale(1.e-6) # convert to µb
    g_pred_fonll_cms_1p45.Scale(1./2.9) # divide by rapidity range

    pt_bins_fonll = np.append(np.asarray(g_pred_fonll_cms_1p45.GetX(), 'd') - np.asarray(g_pred_fonll_cms_1p45.GetEXlow(), 'd'), g_pred_fonll_cms_1p45.GetX()[g_pred_fonll_cms_1p45.GetN()-1] + g_pred_fonll_cms_1p45.GetEXhigh()[g_pred_fonll_cms_1p45.GetN()-1])

    h_pred_fonll_cms_1p45 = ROOT.TH1F('h_pred_fonll', 'h_pred_fonll', g_pred_fonll_cms_1p45.GetN(), pt_bins_fonll)
    for i in range(1, g_pred_fonll_cms_1p45.GetN()+1):
        h_pred_fonll_cms_1p45.SetBinContent(i, g_pred_fonll_cms_1p45.GetY()[i-1])
        h_pred_fonll_cms_1p45.SetBinError(i, 1.e-12)

    # Get predictions
    fonll_file = ROOT.TFile.Open('fonll/fonll_bhadron_nnpdfs_13tev_y2p1.root')
    g_pred_fonll_cms_2p1 = fonll_file.Get('gBhadrNNPDF30')
    fonll_file.Close()

    g_pred_fonll_cms_2p1.Scale(1.e-6) # convert to µb
    g_pred_fonll_cms_2p1.Scale(1./4.2) # divide by rapidity range

    pt_bins_fonll = np.append(np.asarray(g_pred_fonll_cms_2p1.GetX(), 'd') - np.asarray(g_pred_fonll_cms_2p1.GetEXlow(), 'd'), g_pred_fonll_cms_2p1.GetX()[g_pred_fonll_cms_2p1.GetN()-1] + g_pred_fonll_cms_2p1.GetEXhigh()[g_pred_fonll_cms_2p1.GetN()-1])

    h_pred_fonll_cms_2p1 = ROOT.TH1F('h_pred_fonll', 'h_pred_fonll', g_pred_fonll_cms_2p1.GetN(), pt_bins_fonll)
    for i in range(1, g_pred_fonll_cms_2p1.GetN()+1):
        h_pred_fonll_cms_2p1.SetBinContent(i, g_pred_fonll_cms_2p1.GetY()[i-1])
        h_pred_fonll_cms_2p1.SetBinError(i, 1.e-12)

    # Set the style
    colors, _ = root_colors_from_matplotlib_colormap('tab10')
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

    g_pred_fonll_alice.SetLineColorAlpha(colors[0], 0.5)
    g_pred_fonll_alice.SetLineWidth(2)
    g_pred_fonll_alice.SetFillStyle(1001)
    g_pred_fonll_alice.SetFillColorAlpha(colors[0], 0.5)

    h_pred_fonll_alice.SetLineColorAlpha(colors[0], 0.5)
    h_pred_fonll_alice.SetLineWidth(2)
    h_pred_fonll_alice.SetMarkerColorAlpha(colors[0], 0.5)

    g_pred_fonll_cms_1p45.SetLineColorAlpha(colors[1], 0.5)
    g_pred_fonll_cms_1p45.SetLineWidth(2)
    g_pred_fonll_cms_1p45.SetFillStyle(1001)
    g_pred_fonll_cms_1p45.SetFillColorAlpha(colors[1], 0.5)

    h_pred_fonll_cms_1p45.SetLineColorAlpha(colors[1], 0.5)
    h_pred_fonll_cms_1p45.SetLineWidth(2)
    h_pred_fonll_cms_1p45.SetMarkerColorAlpha(colors[1], 0.5)

    g_pred_fonll_cms_2p1.SetLineColorAlpha(colors[3], 0.5)
    g_pred_fonll_cms_2p1.SetLineWidth(2)
    g_pred_fonll_cms_2p1.SetFillStyle(1001)
    g_pred_fonll_cms_2p1.SetFillColorAlpha(colors[3], 0.5)

    h_pred_fonll_cms_2p1.SetLineColorAlpha(colors[3], 0.5)
    h_pred_fonll_cms_2p1.SetLineWidth(2)
    h_pred_fonll_cms_2p1.SetMarkerColorAlpha(colors[3], 0.5)

    # Draw
    c = ROOT.TCanvas('c', 'c', 600, 600)

    pad_axis_label = ROOT.TPad('pad_axis_label', 'pad_axis_label', 0., 0, 1, 1)
    pad_axis_label.SetLogx()
    pad_axis_label.Draw()
    pad_axis_label.cd()
    h_frame_axis_label = pad_axis_label.DrawFrame(1, 0.123, 100, 0.9876, ';#it{p}_{T} (GeV/#it{c});')
    h_frame_axis_label.GetXaxis().SetLabelOffset(0.)
    h_frame_axis_label.GetXaxis().SetTitleOffset(1.1)
    h_frame_axis_label.GetYaxis().SetTitleOffset(1.3)
    h_frame_axis_label.GetXaxis().SetTitleSize(0.035)
    h_frame_axis_label.GetYaxis().SetTitleSize(0.035)
    h_frame_axis_label.GetXaxis().SetLabelSize(0.03)
    h_frame_axis_label.GetYaxis().SetLabelSize(0.03)

    c.cd()
    total_pad = ROOT.TPad('total_pad', 'total_pad', 0., 0.1, 1, 0.95)
    total_pad.Draw()
    total_pad.cd()

    pad_cross_sec = ROOT.TPad('pad_cross_sec', 'pad_cross_sec', 0, 0.2, 1, 1)
    pad_cross_sec.SetBottomMargin(0)
    pad_cross_sec.SetTopMargin(0)
    pad_cross_sec.Draw()
    pad_cross_sec.cd()
    pad_cross_sec.SetLogy()
    pad_cross_sec.SetLogx()
    h_frame = pad_cross_sec.DrawFrame(1, 5e-6, 100, 4.e2, ';#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} (#mub GeV^{#minus1}#kern[0.25]{#it{c}}) ')
    h_frame.GetXaxis().SetTitleOffset(1.1)
    h_frame.GetYaxis().SetTitleOffset(1.195)
    h_frame.GetXaxis().SetTitleSize(0.045)
    h_frame.GetYaxis().SetTitleSize(0.045)
    h_frame.GetXaxis().SetLabelSize(0.045)
    h_frame.GetYaxis().SetLabelSize(0.045)
    g_pred_fonll_cms_2p1.Draw('same 2')
    h_pred_fonll_cms_2p1.Draw('same e')
    g_pred_fonll_alice.Draw('same 2')
    h_pred_fonll_alice.Draw('same e')
    g_pred_fonll_cms_1p45.Draw('same 2')
    h_pred_fonll_cms_1p45.Draw('same e')
    g_stat_cms_1p45.Draw('same pZ')
    g_stat_cms_2p1.Draw('same pZ')
    g_syst_cms_1p45.Draw('same 5')
    g_syst_cms_2p1.Draw('same 5')
    h_stat.Draw('same p')
    g_syst.Draw('same 5')

    # Legend
    leg_alice = ROOT.TLegend(0.6, 0.78, 0.7, 0.88)
    leg_alice.SetBorderSize(0)
    leg_alice.SetFillStyle(0)
    leg_alice.SetTextSize(0.05)
    leg_alice.SetMargin(0.7)
    leg_alice.SetHeader('ALICE, #sqrt{#it{s}} = 13.6 TeV')
    leg_alice.AddEntry(h_stat, '|#it{y}| < 0.5', 'lp')
    leg_alice.Draw()

    leg_cms = ROOT.TLegend(0.14, 0.03, 0.24, 0.18)
    leg_cms.SetBorderSize(0)
    leg_cms.SetFillStyle(0)
    leg_cms.SetTextSize(0.05)
    leg_cms.SetMargin(0.7)
    leg_cms.SetHeader('CMS, #sqrt{#it{s}} = 13 TeV')
    leg_cms.AddEntry(g_stat_cms_1p45, '|#it{y}| < 1.45', 'lp')
    leg_cms.AddEntry(g_stat_cms_2p1, '|#it{y}| < 2.1', 'lp')
    leg_cms.Draw()

    leg_fonll = ROOT.TLegend(0.45, 0.03, 0.55, 0.23)
    leg_fonll.SetBorderSize(0)
    leg_fonll.SetFillStyle(0)
    leg_fonll.SetTextSize(0.05)
    leg_fonll.SetMargin(0.7)
    leg_fonll.SetHeader('FONLL')
    leg_fonll.AddEntry(g_pred_fonll_alice, '#sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5', 'lf')
    leg_fonll.AddEntry(g_pred_fonll_cms_1p45, '#sqrt{#it{s}} = 13 TeV, |#it{y}| < 1.45', 'lf')
    leg_fonll.AddEntry(g_pred_fonll_cms_2p1, '#sqrt{#it{s}} = 13 TeV, |#it{y}| < 2.1', 'lf')
    leg_fonll.Draw()

    # Add the text
    text_decay_alice = ROOT.TLatex(0.61, 0.9, 'B^{0} mesons')
    text_decay_alice.SetNDC()
    text_decay_alice.SetTextSize(0.05)
    text_decay_alice.SetTextFont(42)
    text_decay_alice.Draw()

    text_decay_cms = ROOT.TLatex(0.15, 0.2, 'B^{+} mesons')
    text_decay_cms.SetNDC()
    text_decay_cms.SetTextSize(0.05)
    text_decay_cms.SetTextFont(42)
    text_decay_cms.Draw()

    text_ALICE = ROOT.TLatex(0.15, 0.9, 'ALICE Preliminary')
    text_ALICE.SetNDC()
    text_ALICE.SetTextSize(0.06)
    text_ALICE.SetTextFont(42)
    text_ALICE.Draw()

    text_pp = ROOT.TLatex(0.15, 0.85, 'pp collisions')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.05)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    ROOT.gPad.RedrawAxis()

    # total_pad.cd()
    # pad_ratio_alice = ROOT.TPad('pad_ratio_alice', 'pad_ratio_alice', 0, 0.2, 1, 0.4)
    # pad_ratio_alice.SetLogx()
    # pad_ratio_alice.SetBottomMargin(0)
    # pad_ratio_alice.SetTopMargin(0)
    # pad_ratio_alice.Draw()
    # pad_ratio_alice.cd()
    # h_frame_ratio_alice = pad_ratio_alice.DrawFrame(1, 0.2, 100, 1.8, ';#it{p}_{T} (GeV/#it{c});ALICE/FONLL_{ }')
    # h_frame_ratio_alice.GetXaxis().SetTitleOffset(1.1)
    # h_frame_ratio_alice.GetYaxis().SetTitleOffset(0.4)
    # h_frame_ratio_alice.GetXaxis().SetTitleSize(0.12)
    # h_frame_ratio_alice.GetYaxis().SetTitleSize(0.13)
    # h_frame_ratio_alice.GetXaxis().SetLabelSize(0.12)
    # h_frame_ratio_alice.GetYaxis().SetLabelSize(0.12)
    # h_frame_ratio_alice.GetYaxis().SetNdivisions(505)

    total_pad.cd()
    pad_ratio_cms = ROOT.TPad('pad_ratio_cms', 'pad_ratio_cms', 0, 0., 1, 0.2)
    pad_ratio_cms.SetBottomMargin(0)
    pad_ratio_cms.SetTopMargin(0)
    pad_ratio_cms.SetLogx()
    pad_ratio_cms.Draw()
    pad_ratio_cms.cd()
    h_frame_ratio_cms = pad_ratio_cms.DrawFrame(1, 0.2, 100, 2.75, ';#it{p}_{T} (GeV/#it{c});Data/FONLL_{ }')
    h_frame_ratio_cms.GetXaxis().SetTitleOffset(1.1)
    h_frame_ratio_cms.GetYaxis().SetTitleOffset(0.33)
    h_frame_ratio_cms.GetYaxis().SetLabelOffset(0.01)
    h_frame_ratio_cms.GetXaxis().SetTitleSize(0.17)
    h_frame_ratio_cms.GetYaxis().SetTitleSize(0.17)
    h_frame_ratio_cms.GetXaxis().SetLabelSize(0.17)
    h_frame_ratio_cms.GetYaxis().SetLabelSize(0.17)
    h_frame_ratio_cms.GetYaxis().SetNdivisions(505)
    h_frame_ratio_cms.GetXaxis().SetTickLength(0.12)

    g_fonll_unc_alice = g_pred_fonll_alice.Clone('g_fonll_unc_alice')
    for i in range(0, g_fonll_unc_alice.GetN()):
        g_fonll_unc_alice.SetPoint(i, g_pred_fonll_alice.GetX()[i], 1)
        g_fonll_unc_alice.SetPointError(i, g_pred_fonll_alice.GetErrorXlow(i), g_pred_fonll_alice.GetErrorXhigh(i), g_pred_fonll_alice.GetErrorYlow(i)/g_pred_fonll_alice.GetY()[i], g_pred_fonll_alice.GetErrorYhigh(i)/g_pred_fonll_alice.GetY()[i])

    g_fonll_unc_alice.Draw("same 2")

    line_one_alice = ROOT.TLine(1, 1, 100, 1)
    line_one_alice.SetLineStyle(2)
    line_one_alice.SetLineColor(ROOT.kBlack)
    line_one_alice.Draw("same")

    h_ratio_data_fonll_stat = h_stat.Clone('h_ratio_data_fonll_stat')
    h_ratio_data_fonll_stat.Divide(h_pred_fonll_alice)

    g_ratio_data_fonll_syst = g_syst.Clone('g_ratio_data_fonll_syst')
    for i in range(0, g_ratio_data_fonll_syst.GetN()):
        g_ratio_data_fonll_syst.SetPoint(i, g_syst.GetX()[i], g_syst.GetY()[i]/g_pred_fonll_alice.GetY()[i])
        g_ratio_data_fonll_syst.SetPointError(i, g_syst.GetErrorX(i), g_syst.GetErrorY(i)/g_pred_fonll_alice.GetY()[i])



    g_fonll_unc_cms_1p45 = g_pred_fonll_cms_1p45.Clone('g_fonll_unc_cms_1p45')
    for i in range(0, g_fonll_unc_cms_1p45.GetN()):
        g_fonll_unc_cms_1p45.SetPoint(i, g_pred_fonll_cms_1p45.GetX()[i], 1)
        g_fonll_unc_cms_1p45.SetPointError(i, g_pred_fonll_cms_1p45.GetErrorXlow(i), g_pred_fonll_cms_1p45.GetErrorXhigh(i), g_pred_fonll_cms_1p45.GetErrorYlow(i)/g_pred_fonll_cms_1p45.GetY()[i], g_pred_fonll_cms_1p45.GetErrorYhigh(i)/g_pred_fonll_cms_1p45.GetY()[i])

    g_fonll_unc_cms_1p45.Draw("same 2")

    line_one_cms = ROOT.TLine(1, 1, 100, 1)
    line_one_cms.SetLineStyle(2)
    line_one_cms.SetLineColor(ROOT.kBlack)
    line_one_cms.Draw("same")

    g_ratio_data_fonll_stat_1p45 = g_stat_cms_1p45.Clone('g_ratio_data_fonll_stat_1p45')
    for i in range(0, g_ratio_data_fonll_stat_1p45.GetN()):
        g_ratio_data_fonll_stat_1p45.SetPoint(i, g_stat_cms_1p45.GetX()[i], g_stat_cms_1p45.GetY()[i]/g_pred_fonll_cms_1p45.GetY()[i])
        g_ratio_data_fonll_stat_1p45.SetPointError(i, g_stat_cms_1p45.GetErrorXlow(i), g_stat_cms_1p45.GetErrorXhigh(i), g_stat_cms_1p45.GetErrorYlow(i)/g_pred_fonll_cms_1p45.GetY()[i], g_stat_cms_1p45.GetErrorYhigh(i)/g_pred_fonll_cms_1p45.GetY()[i])
    g_ratio_data_fonll_stat_1p45.Draw('same pZ')

    g_ratio_data_fonll_syst_1p45 = g_syst_cms_1p45.Clone('g_ratio_data_fonll_syst_1p45')
    for i in range(0, g_ratio_data_fonll_syst_1p45.GetN()):
        g_ratio_data_fonll_syst_1p45.SetPoint(i, g_syst_cms_1p45.GetX()[i], g_syst_cms_1p45.GetY()[i]/g_pred_fonll_cms_1p45.GetY()[i])
        g_ratio_data_fonll_syst_1p45.SetPointError(i, g_syst_cms_1p45.GetErrorXlow(i), g_syst_cms_1p45.GetErrorXhigh(i), g_syst_cms_1p45.GetErrorYlow(i)/g_pred_fonll_cms_1p45.GetY()[i], g_syst_cms_1p45.GetErrorYhigh(i)/g_pred_fonll_cms_1p45.GetY()[i])

    g_ratio_data_fonll_syst_1p45.Draw('same 5') 

    g_fonll_unc_cms_2p1 = g_pred_fonll_cms_2p1.Clone('g_fonll_unc_cms_2p1')
    for i in range(0, g_fonll_unc_cms_2p1.GetN()):
        g_fonll_unc_cms_2p1.SetPoint(i, g_pred_fonll_cms_2p1.GetX()[i], 1)
        g_fonll_unc_cms_2p1.SetPointError(i, g_pred_fonll_cms_2p1.GetErrorXlow(i), g_pred_fonll_cms_2p1.GetErrorXhigh(i), g_pred_fonll_cms_2p1.GetErrorYlow(i)/g_pred_fonll_cms_2p1.GetY()[i], g_pred_fonll_cms_2p1.GetErrorYhigh(i)/g_pred_fonll_cms_2p1.GetY()[i])

    g_fonll_unc_cms_2p1.Draw("same 2")

    g_ratio_data_fonll_stat_2p1 = g_stat_cms_2p1.Clone('g_ratio_data_fonll_stat_2p1')
    for i in range(0, g_ratio_data_fonll_stat_2p1.GetN()):
        g_ratio_data_fonll_stat_2p1.SetPoint(i, g_stat_cms_2p1.GetX()[i], g_stat_cms_2p1.GetY()[i]/g_pred_fonll_cms_2p1.GetY()[i])
        g_ratio_data_fonll_stat_2p1.SetPointError(i, g_stat_cms_2p1.GetErrorXlow(i), g_stat_cms_2p1.GetErrorXhigh(i), g_stat_cms_2p1.GetErrorYlow(i)/g_pred_fonll_cms_2p1.GetY()[i], g_stat_cms_2p1.GetErrorYhigh(i)/g_pred_fonll_cms_2p1.GetY()[i])
    g_ratio_data_fonll_stat_2p1.Draw('same pZ')

    g_ratio_data_fonll_syst_2p1 = g_syst_cms_2p1.Clone('g_ratio_data_fonll_syst_2p1')
    for i in range(0, g_ratio_data_fonll_syst_2p1.GetN()):
        g_ratio_data_fonll_syst_2p1.SetPoint(i, g_syst_cms_2p1.GetX()[i], g_syst_cms_2p1.GetY()[i]/g_pred_fonll_cms_2p1.GetY()[i])
        g_ratio_data_fonll_syst_2p1.SetPointError(i, g_syst_cms_2p1.GetErrorXlow(i), g_syst_cms_2p1.GetErrorXhigh(i), g_syst_cms_2p1.GetErrorYlow(i)/g_pred_fonll_cms_2p1.GetY()[i], g_syst_cms_2p1.GetErrorYhigh(i)/g_pred_fonll_cms_2p1.GetY()[i])

    g_ratio_data_fonll_syst_2p1.Draw('same 5') 

    h_ratio_data_fonll_stat.Draw('same p')
    g_ratio_data_fonll_syst.Draw('same 5') 

    ROOT.gPad.RedrawAxis()

    c.SaveAs('figures/cross_section/cross_section_vs_CMS.pdf')

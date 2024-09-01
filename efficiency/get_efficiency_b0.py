'''
Script for the computation of the B0 meson efficiency
run: python get_efficiency_b0.py config_file_name.yml
'''

import ROOT 
import pandas as pd 
import numpy as np
import argparse
import yaml
import ctypes
import sys
sys.path.append('Utils')
from DfUtils import read_parquet_in_batches
from AnalysisUtils import evaluate_efficiency_from_histos

def draw_efficiency_figure(h_eff, h_acc, out_file_name_pdf):
    """
    Draw the efficiency and acceptance histograms.

    Args:
        h_eff (TH1F): The efficiency histogram.
        h_acc (TH1F): The acceptance histogram.
        out_file_name_pdf (str): The name of the output PDF file.

    Returns:
        None
    """

    # Set style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    h_eff.SetMarkerStyle(ROOT.kFullCircle)
    h_eff.SetMarkerColor(ROOT.kAzure+2)
    h_eff.SetMarkerSize(2)
    h_eff.SetLineColor(ROOT.kAzure+2)
    h_eff.SetLineWidth(2)

    h_acc.SetMarkerStyle(ROOT.kFullDiamond)
    h_acc.SetMarkerColor(ROOT.kTeal-7)
    h_acc.SetMarkerSize(2.5)
    h_acc.SetLineColor(ROOT.kTeal-7)
    h_acc.SetLineWidth(2)
    h_acc.SetLineStyle(7)

    leg = ROOT.TLegend(0.5, 0.2, 0.8, 0.3)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(h_acc, "Acceptance", "pl")
    leg.AddEntry(h_eff, "Acc. #times #varepsilon", "pl")

    c_eff = ROOT.TCanvas('c_eff', '', 800, 800)
    c_eff.SetLogy()
    h_frame = c_eff.DrawFrame(h_eff.GetBinLowEdge(1), 1.e-4, h_eff.GetBinLowEdge(h_eff.GetNbinsX()+1), 10.,
                ';#it{p}_{T} (GeV/#it{c});Acceptance #times Efficiency;')
    h_frame.GetXaxis().SetTitleOffset(1.2)
    h_frame.GetYaxis().SetTitleOffset(1.5)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_eff.Draw('][ hist same')
    h_eff.DrawClone('pe X0 same')
    h_acc.Draw('][ hist same')
    h_acc.DrawClone('pe X0 same')
    leg.Draw()

    # Add the text
    text_decay = ROOT.TLatex(0.51, 0.36, 'B^{0}#rightarrow D^{#font[122]{-}}#pi^{+}#rightarrow #pi^{#font[122]{-}}K^{+}#pi^{#font[122]{-}}#pi^{+}')
    text_decay.SetNDC()
    text_decay.SetTextSize(0.04)
    text_decay.SetTextFont(42)
    text_decay.Draw()

    text_conj = ROOT.TLatex(0.51, 0.31, 'and charge conjugate')
    text_conj.SetNDC()
    text_conj.SetTextSize(0.04)
    text_conj.SetTextFont(42)
    text_conj.Draw()

    text_alice = ROOT.TLatex(0.18, 0.88, 'Work in Progress')
    text_alice.SetNDC()
    text_alice.SetTextSize(0.06)
    text_alice.SetTextFont(42)
    text_alice.Draw()

    text_pp = ROOT.TLatex(0.18, 0.83, 'pp collisions, #sqrt{#it{s}} = 13.6 TeV')
    text_pp.SetNDC()
    text_pp.SetTextSize(0.04)
    text_pp.SetTextFont(42)
    text_pp.Draw()

    text_rapidity = ROOT.TLatex(0.18, 0.78, '|y| < 0.5')
    text_rapidity.SetNDC()
    text_rapidity.SetTextSize(0.04)
    text_rapidity.SetTextFont(42)
    text_rapidity.Draw()

    c_eff.SaveAs(out_file_name_pdf)

def compute_efficiency(config_file_name):
    """
    Compute the efficiency of a B0 particle based on the given configuration and cut set.

    Args:
        config_file_name (str): The file name of the configuration file.

    Returns:
        None
    """

    with open(config_file_name, 'r') as yml_config_file:
        config = yaml.safe_load(yml_config_file)
    
    with open(config['cutset_file_name'], 'r') as yml_cut_set_file:
        cut_set = yaml.safe_load(yml_cut_set_file)

    pt_mins = cut_set['pt']['mins']
    pt_maxs = cut_set['pt']['maxs']
    pt_lims = list(pt_mins)
    n_pt_bins = len(pt_mins)
    pt_lims.append(pt_maxs[-1])

    h_reco_integrated = ROOT.TH1F('h_reco_integrated', ';#it{p}_{T} (GeV/#it{c});Reconstructed', n_pt_bins, np.asarray(pt_lims, 'd'))
    h_gen_integrated = ROOT.TH1F('h_gen_integrated', ';#it{p}_{T} (GeV/#it{c});Generated', n_pt_bins, np.asarray(pt_lims, 'd'))
    h_gen_in_acc_integrated = ROOT.TH1F('h_gen_in_acc_integrated', ';#it{p}_{T} (GeV/#it{c});Generated in acceptance', n_pt_bins, np.asarray(pt_lims, 'd'))

    # Get the generated and reconstructed particles
    for i_file, in_file_name in enumerate(config['gen']['file_names']):
        in_file = ROOT.TFile.Open(in_file_name)
        h_sparse_gen = in_file.Get(config['gen']['sparse_name'])
        h_sparse_gen_in_acc = in_file.Get(config['gen']['sparse_name_acc'])
        if i_file == 0:
            h_gen = h_sparse_gen.Projection(config['gen']['pt_axis'])
            h_gen.SetDirectory(0)
            h_gen_in_acc = h_sparse_gen_in_acc.Projection(config['gen']['pt_axis'])
            h_gen_in_acc.SetDirectory(0)
        else:
            h_gen.Add(h_sparse_gen.Projection(config['gen']['pt_axis']))
            h_gen_in_acc.Add(h_sparse_gen_in_acc.Projection(config['gen']['pt_axis']))
        in_file.Close()

    h_gen.SetName('h_gen')
    h_gen_in_acc.SetName('h_gen_in_acc')
    h_reco = h_gen.Clone('h_reco')
    h_reco.Reset()


    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):

        # Apply the cuts on the reconstructed particles
        df_reco = pd.concat([read_parquet_in_batches(parquet, f"{pt_min} < fPt < {pt_max}") for parquet in config['reco_file_names']])
        
        sel_to_apply = ''
        for cut_var in cut_set:
            if cut_var == 'pt' or cut_var == 'M':
                continue
            sel_to_apply += f" {cut_set[cut_var]['mins'][i_pt]} < {cut_var} < {cut_set[cut_var]['maxs'][i_pt]} and"
        sel_to_apply = sel_to_apply[:-3] # Remove the last 'and'

        df_reco = df_reco.query(sel_to_apply)
        for pt in df_reco['fPt']:
            h_reco.Fill(pt)

        # Get the number of generated and reconstructed particles in the given pt range
        n_reco_unc, n_gen_unc, n_gen_in_acc_unc = (ctypes.c_double() for _ in range(3))
        n_reco = h_reco.IntegralAndError(h_reco.FindBin(pt_min), h_reco.FindBin(pt_max)-1, n_reco_unc)
        n_gen = h_gen.IntegralAndError(h_gen.FindBin(pt_min), h_gen.FindBin(pt_max)-1, n_gen_unc)
        n_gen_in_acc = h_gen_in_acc.IntegralAndError(h_gen_in_acc.FindBin(pt_min), h_gen_in_acc.FindBin(pt_max)-1, n_gen_in_acc_unc)

        h_reco_integrated.SetBinContent(i_pt+1, n_reco)
        h_reco_integrated.SetBinError(i_pt+1, n_reco_unc.value)
        h_gen_integrated.SetBinContent(i_pt+1, n_gen)
        h_gen_integrated.SetBinError(i_pt+1, n_gen_unc.value)
        h_gen_in_acc_integrated.SetBinContent(i_pt+1, n_gen_in_acc)
        h_gen_in_acc_integrated.SetBinError(i_pt+1, n_gen_in_acc_unc.value)

    # Compute the efficiency and acceptance in the given pt ranges
    h_eff = evaluate_efficiency_from_histos(h_gen_integrated, h_reco_integrated)
    h_eff.SetName('h_eff')

    h_acc = evaluate_efficiency_from_histos(h_gen_integrated, h_gen_in_acc_integrated)
    h_acc.SetName('h_acc')

    # Compute the efficiency in pt bins used for the generated particles
    h_eff_fine_bins = evaluate_efficiency_from_histos(h_gen, h_reco)
    h_eff_fine_bins.SetName('h_eff_fine_bins')

    h_acc_fine_bins = evaluate_efficiency_from_histos(h_gen, h_gen_in_acc)
    h_acc_fine_bins.SetName('h_acc_fine_bins')

    out_file = ROOT.TFile(config['output_file_name'], 'recreate')
    h_eff.Write()
    h_eff_fine_bins.Write()
    h_acc.Write()
    h_acc_fine_bins.Write()
    h_gen.Write()
    h_reco.Write()
    h_gen_integrated.Write()
    h_reco_integrated.Write()
    out_file.Close()

    out_file_name_pdf = config['output_file_name'].replace('.root', '.pdf')
    draw_efficiency_figure(h_eff, h_acc, out_file_name_pdf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('config_file_name', metavar='text', default='config_eff.yaml')
    args = parser.parse_args()

    compute_efficiency(args.config_file_name)

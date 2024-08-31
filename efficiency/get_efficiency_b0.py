'''
Script for the computation of the B0 meson efficiency
run: python get_efficiency_b0.py config_file_name.yml cut_set_file_name.yml
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

def compute_efficiency(config_file_name, cut_set_file_name):
    """
    Compute the efficiency of a B0 particle based on the given configuration and cut set.

    Args:
        config_file_name (str): The file name of the configuration file.
        cut_set_file_name (str): The file name of the cut set file.

    Returns:
        None
    """

    with open(args.config_file_name, 'r') as yml_config_file:
        config = yaml.safe_load(yml_config_file)
    
    with open(config['cutset_file_name'], 'r') as yml_cut_set_file:
        cut_set = yaml.safe_load(yml_cut_set_file)

    # Set style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetTitleOffset(1.2, 'X')
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    pt_mins = cut_set['pt']['mins']
    pt_maxs = cut_set['pt']['maxs']
    pt_lims = list(pt_mins)
    n_pt_bins = len(pt_mins)
    pt_lims.append(pt_maxs[-1])

    h_reco_integrated = ROOT.TH1F('h_reco_integrated', ';#it{p}_{T} (GeV/#it{c});Reconstructed', n_pt_bins, np.asarray(pt_lims, 'd'))
    h_gen_integrated = ROOT.TH1F('h_gen_integrated', ';#it{p}_{T} (GeV/#it{c});Generated', n_pt_bins, np.asarray(pt_lims, 'd'))

    for i_file, in_file_name in enumerate(config['gen']['file_names']):
        in_file = ROOT.TFile.Open(in_file_name)
        h_sparse_gen = in_file.Get(config['gen']['sparse_name'])
        if i_file == 0:
            h_gen = h_sparse_gen.Projection(config['gen']['pt_axis'])
            h_gen.SetDirectory(0)
        else:
            h_gen.Add(h_sparse_gen.Projection(config['gen']['pt_axis']))
        in_file.Close()

    h_gen.SetName('h_gen')
    h_reco = h_gen.Clone('h_reco')
    h_reco.Reset()


    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):

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
        n_reco_unc, n_gen_unc = (ctypes.c_double() for _ in range(2))
        n_reco = h_reco.IntegralAndError(h_reco.FindBin(pt_min), h_reco.FindBin(pt_max)-1, n_reco_unc)
        n_gen = h_gen.IntegralAndError(h_gen.FindBin(pt_min), h_gen.FindBin(pt_max)-1, n_gen_unc)

        h_reco_integrated.SetBinContent(i_pt+1, n_reco)
        h_reco_integrated.SetBinError(i_pt+1, n_reco_unc.value)
        h_gen_integrated.SetBinContent(i_pt+1, n_gen)
        h_gen_integrated.SetBinError(i_pt+1, n_gen_unc.value)

    # Compute the efficiency in the given pt ranges
    h_eff = evaluate_efficiency_from_histos(h_gen_integrated, h_reco_integrated)
    h_eff.SetMarkerStyle(ROOT.kFullCircle)
    h_eff.SetMarkerColor(ROOT.kAzure+3)
    h_eff.SetMarkerSize(1.5)
    h_eff.SetLineColor(ROOT.kAzure+3)
    h_eff.SetLineWidth(2)

    # Compute the efficiency in pt bins used for the generated particles
    h_eff_fine_bins = evaluate_efficiency_from_histos(h_gen, h_reco)
    h_eff_fine_bins.SetName('h_eff_fine_bins')

    leg = ROOT.TLegend(0.6, 0.2, 0.8, 0.4)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(h_eff, " B^{0}", "p")

    c_eff = ROOT.TCanvas('c_eff', '', 800, 800)
    c_eff.DrawFrame(0, 1.e-4, pt_maxs[n_pt_bins-1], 1.,
                ';#it{p}_{T} (GeV/#it{c});Efficiency;')
    c_eff.SetLogy()
    h_eff.Draw('same')
    leg.Draw()

    out_file = ROOT.TFile(config['output_file_name'], 'recreate')
    h_eff.Write()
    h_eff_fine_bins.Write()
    h_gen.Write()
    h_reco.Write()
    h_gen_integrated.Write()
    h_reco_integrated.Write()
    out_file.Close()

    out_file_name_pdf = config['output_file_name'].replace('.root', '.pdf')
    c_eff.SaveAs(out_file_name_pdf)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('config_file_name', metavar='text', default='config_Efficiency')
    args = parser.parse_args()

    compute_efficiency(args.config_file_name)

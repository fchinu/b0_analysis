import argparse
import os

import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt

import uproot
import ROOT

def get_integrated_cross_section(config_filename):
    with open(config_filename) as f:
        config = yaml.safe_load(f)

    # Get the measured cross section
    correlated_uncs = {}
    uncorrelated_uncs = {}
    with uproot.open(config['inputs']['cross_section']) as f:
        h_stat = f['h_stat']
        for syst in config['systematics']['correlated']:
            correlated_uncs[syst] = f[f'{syst}']
        for syst in config['systematics']['uncorrelated']:
            uncorrelated_uncs[syst] = f[f'{syst}']
    pt_bins = h_stat.axis().edges()
    bin_widths = np.diff(pt_bins)

    # Get the visible cross section
    visible_cross_section = np.sum(h_stat.values()*bin_widths)
    visible_cross_section_unc_stat = np.sqrt(np.sum(np.square(h_stat.errors() * bin_widths)))
    visible_cross_section_unc_syst = 0
    corr_syst_unc = {}
    uncorr_syst_unc = {}
    for syst in correlated_uncs:
        corr_syst_unc[syst] = np.sum(np.abs(correlated_uncs[syst].values() * bin_widths))
        visible_cross_section_unc_syst += np.sum(np.abs(correlated_uncs[syst].values() * bin_widths)) ** 2
    for syst in uncorrelated_uncs:
        uncorr_syst_unc[syst] = np.sqrt(np.sum(np.square(uncorrelated_uncs[syst].values() * bin_widths)))
        visible_cross_section_unc_syst += np.sum(np.square(uncorrelated_uncs[syst].values() * bin_widths))

    h_stat_vis = ROOT.TH1D('h_stat_vis', 'h_stat_vis', 1, pt_bins[0], pt_bins[-1])
    h_syst_vis_tot = ROOT.TH1D('h_syst_vis_tot', 'h_syst_vis_tot', 1, pt_bins[0], pt_bins[-1])
    correlated_vis = ROOT.TH1D('correlated_vis', 'correlated_vis', 1, pt_bins[0], pt_bins[-1])
    uncorrelated_vis = ROOT.TH1D('uncorrelated_vis', 'uncorrelated_vis', 1, pt_bins[0], pt_bins[-1])
    h_systs_vis = {}
    h_systs_vis_rel = {}
    for syst in correlated_uncs:
        h_systs_vis[syst] = ROOT.TH1D(syst, syst, 1, pt_bins[0], pt_bins[-1])
        h_systs_vis_rel[syst] = ROOT.TH1D(f'{syst}_rel', f'{syst}_rel', 1, pt_bins[0], pt_bins[-1])
    for syst in uncorrelated_uncs:
        h_systs_vis[syst] = ROOT.TH1D(syst, syst, 1, pt_bins[0], pt_bins[-1])
        h_systs_vis_rel[syst] = ROOT.TH1D(f'{syst}_rel', f'{syst}_rel', 1, pt_bins[0], pt_bins[-1])

    h_stat_vis.SetBinContent(1, visible_cross_section)
    h_stat_vis.SetBinError(1, visible_cross_section_unc_stat)

    h_syst_vis_tot.SetBinContent(1, visible_cross_section)
    h_syst_vis_tot.SetBinError(1, np.sqrt(visible_cross_section_unc_syst))

    correlated_vis.SetBinContent(1, np.sqrt(np.sum(np.square(list(corr_syst_unc.values())))))
    uncorrelated_vis.SetBinContent(1, np.sqrt(np.sum(np.square(list(uncorr_syst_unc.values())))))
    for syst in corr_syst_unc:
        h_systs_vis[syst].SetBinContent(1, corr_syst_unc[syst])
        h_systs_vis_rel[syst].SetBinContent(1, corr_syst_unc[syst]/visible_cross_section)
    for syst in uncorr_syst_unc:
        h_systs_vis[syst].SetBinContent(1, uncorr_syst_unc[syst])
        h_systs_vis_rel[syst].SetBinContent(1, uncorr_syst_unc[syst]/visible_cross_section) 


    col_names = ["ptmin", "ptmax", "central", "min", "max", "min_sc", "max_sc",
                "min_mass", "max_mass", "min_pdf", "max_pdf", "fr_dot5_dot5", "fr_2_2",
                "fr_2_1", "fr_1_2", "fr_1_dot5", "fr_dot5_1"]
    
    df_fonll_integrated = pd.read_csv(config['inputs']['fonll_integrated'], names=col_names, sep='\s+', comment='#')
    df_fonll_visible = pd.read_csv(config['inputs']['fonll_visible'], names=col_names, sep='\s+', comment='#')
    
    cols_scale_variations = ["central", "min_mass", "max_mass", "min_pdf", "max_pdf", "fr_dot5_dot5", "fr_2_2", "fr_2_1", "fr_1_2", "fr_1_dot5", "fr_dot5_1"]

    # Compute scale factors
    scale_factors = {}
    for scale in cols_scale_variations:
        scale_factors[scale] = float(df_fonll_integrated[scale] / df_fonll_visible[scale])

    # Plotting with markers
    plt.figure(figsize=(10, 6))
    variations = list(scale_factors.keys())
    factors = list(scale_factors.values())

    # Create the scatter plot
    plt.scatter(factors, variations, color='C0', marker='o', label='Scale Factor')

    # Annotate each point with its value
    for factor, variation in zip(factors, variations):
        plt.text(factor+5.e-5, variation, f'{factor:.3f}', va='center', ha='left', fontsize=8)

    # Add labels and title
    plt.xlabel('Scale Factor')
    plt.ylabel('Variation')
    plt.title('Scale Factors for Different Variations')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Save and show the plot
    output_path = os.path.join(config['output']['dir'], 'scale_factors.png')
    plt.savefig(output_path)
    plt.show()

    h_stat_int = ROOT.TH1D('h_stat_int', 'h_stat_int', 1, 0, 1)
    h_syst_int_tot = ROOT.TH1D('h_syst_int_tot', 'h_syst_int_tot', 1, 0, 1)
    correlated_int = ROOT.TH1D('correlated_int', 'correlated_int', 1, 0, 1)
    uncorrelated_int = ROOT.TH1D('uncorrelated_int', 'uncorrelated_int', 1, 0, 1)

    # Set bin content for the new histograms
    h_stat_int.SetBinContent(1, h_stat_vis.GetBinContent(1) * scale_factors['central'])
    h_stat_int.SetBinError(1, h_stat_vis.GetBinError(1) * scale_factors['central'])

    h_syst_int_tot.SetBinContent(1, h_syst_vis_tot.GetBinContent(1) * scale_factors['central'])
    h_syst_int_tot.SetBinError(1, h_syst_vis_tot.GetBinError(1) * scale_factors['central'])

    correlated_int.SetBinContent(1, correlated_vis.GetBinContent(1) * scale_factors['central'])
    correlated_int.SetBinError(1, correlated_vis.GetBinError(1) * scale_factors['central'])

    uncorrelated_int.SetBinContent(1, uncorrelated_vis.GetBinContent(1) * scale_factors['central'])
    uncorrelated_int.SetBinError(1, uncorrelated_vis.GetBinError(1) * scale_factors['central'])

    # Create dictionaries for systematic uncertainties
    h_systs_int = {}
    h_systs_int_rel = {}

    # Loop over correlated uncertainties
    for syst in correlated_uncs:
        h_systs_int[syst] = ROOT.TH1D(syst, syst, 1, 0, 1)
        h_systs_int[syst].SetBinContent(1, h_systs_vis[syst].GetBinContent(1) * scale_factors['central'])
        h_systs_int[syst].SetBinError(1, h_systs_vis[syst].GetBinError(1) * scale_factors['central'])

        h_systs_int_rel[syst] = ROOT.TH1D(f'{syst}_rel', f'{syst}_rel', 1, 0, 1)
        h_systs_int_rel[syst].SetBinContent(1, h_systs_int[syst].GetBinContent(1)/h_stat_int.GetBinContent(1))

    # Loop over uncorrelated uncertainties
    for syst in correlated_uncs:
        h_systs_int[syst] = ROOT.TH1D(syst, syst, 1, 0, 1)
        h_systs_int[syst].SetBinContent(1, h_systs_vis[syst].GetBinContent(1) * scale_factors['central'])
        h_systs_int[syst].SetBinError(1, h_systs_vis[syst].GetBinError(1) * scale_factors['central'])

        h_systs_int_rel[syst] = ROOT.TH1D(f'{syst}_rel', f'{syst}_rel', 1, 0, 1)
        h_systs_int_rel[syst].SetBinContent(1, h_systs_int[syst].GetBinContent(1)/h_stat_int.GetBinContent(1))

    max_extrap_factor = max(scale_factors.values())
    min_extrap_factor = min(scale_factors.values())

    g_extrap = ROOT.TGraphAsymmErrors(h_stat_int)
    g_extrap.SetPoint(0, 0.5, h_stat_int.GetBinContent(1))
    g_extrap.SetPointError(
        0, 0.5, 0.5,
        (scale_factors['central'] - min_extrap_factor) * h_stat_int.GetBinContent(1), 
        (max_extrap_factor - scale_factors['central']) * h_stat_int.GetBinContent(1), 
    )
    
    g_extrap_rel = ROOT.TGraphAsymmErrors(g_extrap)
    g_extrap_rel.Scale(1/h_stat_int.GetBinContent(1))

    out_file = os.path.join(
        config['output']['dir'],
        config['output']['file_name']
    )
    with ROOT.TFile.Open(out_file, 'recreate') as f:
        f.mkdir('visible')
        f.cd('visible')
        h_stat_vis.Write()
        h_syst_vis_tot.Write()
        correlated_vis.Write()
        uncorrelated_vis.Write()
        for syst in h_systs_vis:
            h_systs_vis[syst].Write()
        for syst in h_systs_vis:
            h_systs_vis_rel[syst].Write()

        f.mkdir('integrated')
        f.cd('integrated')
        h_stat_int.Write()
        h_syst_int_tot.Write()
        correlated_int.Write()
        uncorrelated_int.Write()
        for syst in h_systs_int:
            h_systs_int[syst].Write()
        for syst in h_systs_int_rel:
            h_systs_int_rel[syst].Write()
        g_extrap.Write('g_extrap')
        g_extrap_rel.Write('g_extrap_rel')
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to compute the integrated cross section')
    parser.add_argument('config_file', help='Path to the config file')
    args = parser.parse_args()

    get_integrated_cross_section(args.config_file)
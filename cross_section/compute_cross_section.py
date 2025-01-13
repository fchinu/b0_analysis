"""
Script for evaluating the cross section.

Run:
    python compute_cross_section.py config.yaml
"""
import argparse
import sys
sys.path.append('utils') # pylint: disable=wrong-import-position
import ROOT # pylint: disable=import-error
import yaml
from analysis_utils import get_n_events_from_zorro # pylint: disable=import-error
# pylint: disable=no-member

def main(config_file_name):
    '''
    Compute the cross section using the provided configuration file.

    Parameters:
    - config_file_name (str): The path to the configuration file.

    Returns:
    None
    '''

    with open(config_file_name, 'r', encoding='utf-8') as yml_config_file:
        config = yaml.load(yml_config_file, yaml.FullLoader)

    analysis_results_files = config['lumi']['analysis_results_files']
    n_events = get_n_events_from_zorro(
        analysis_results_files, config['lumi']['zorro_folder'],
        config['lumi']['triggers_of_interest'], config['lumi']['h_collisions_path']
    )
    int_lumi_before_bc = n_events/config['lumi']['tvx_cross_section']

    lumi_before_bc, lumi_after_bc = 0, 0
    for file in analysis_results_files:
        infile = ROOT.TFile.Open(file)
        lumi_before_bc += infile.Get(f"{config['lumi']['folder_bc_cuts']}/hLumiTVX").Integral()
        lumi_after_bc += infile.Get(f"{config['lumi']['folder_bc_cuts']}/hLumiTVXafterBCcuts").Integral()
        infile.Close()
    
    int_lumi_after_bc = int_lumi_before_bc * lumi_after_bc / lumi_before_bc

    br_b_to_d = config['br']['b0_todminuspi']
    br_d_to_pikpi = config['br']['dplus_tokpipi']

    # Get the raw yields
    raw_yield_file = ROOT.TFile.Open(config['rawyield_file'])
    h_raw_yield = raw_yield_file.Get('h_rawyields')
    h_raw_yield.SetDirectory(0)
    raw_yield_file.Close()

    # Get the efficiency
    eff_file = ROOT.TFile.Open(config['efficiency_file'])
    h_eff = eff_file.Get('h_eff')
    h_eff.SetDirectory(0)
    eff_file.Close()

    # Compute the cross section
    h_cross_section = h_raw_yield.Clone('h_cross_section')
    h_cross_section.Divide(h_eff)
    h_cross_section.Scale(1/(2 * br_b_to_d * br_d_to_pikpi * int_lumi_after_bc))
    for i in range(1, h_cross_section.GetNbinsX()+1):
        h_cross_section.SetBinContent(i, h_cross_section.GetBinContent(i) / h_cross_section.GetBinWidth(i))
        h_cross_section.SetBinError(i, h_cross_section.GetBinError(i) / h_cross_section.GetBinWidth(i))

    h_lumi_before_bc = h_cross_section.Clone()
    h_lumi_after_bc = h_cross_section.Clone()
    for i in range(1, h_cross_section.GetNbinsX()+1):
        h_lumi_before_bc.SetBinContent(i, int_lumi_before_bc)
        h_lumi_before_bc.SetBinError(i, 0.1 * int_lumi_before_bc)
        h_lumi_after_bc.SetBinContent(i, int_lumi_after_bc)
        h_lumi_after_bc.SetBinError(i, 0.1 * int_lumi_after_bc)

    outfile = ROOT.TFile.Open(config['output_file'], 'recreate')
    h_cross_section.Write()
    h_lumi_before_bc.Write('h_lumi_before_bc')
    h_lumi_after_bc.Write('h_lumi_after_bc')
    outfile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('configFileName', metavar='text', default='config.yaml')
    args = parser.parse_args()
    main(args.configFileName)

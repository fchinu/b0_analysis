import argparse 
import ROOT
import yaml
import sys
sys.path.append('Utils')
from AnalysisUtils import get_n_events_from_zorro

def main(config_file_name):
    '''
    Compute the cross section using the provided configuration file.

    Parameters:
    - config_file_name (str): The path to the configuration file.

    Returns:
    None
    '''
    
    with open(config_file_name, 'r') as yml_config_file:
        config = yaml.load(yml_config_file, yaml.FullLoader)

    analysis_results_files = config['lumi']['analysis_results_files']
    n_events = get_n_events_from_zorro(analysis_results_files, config['lumi']['zorro_folder'], config['lumi']['triggers_of_interest'], config['lumi']['h_collisions_path'])
    int_lumi = n_events/config['lumi']['tvx_cross_section']

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
    h_cross_section.Scale(1/(2 * br_b_to_d * br_d_to_pikpi * int_lumi))
    for i in range(1, h_cross_section.GetNbinsX()+1):
        h_cross_section.SetBinContent(i, h_cross_section.GetBinContent(i) / h_cross_section.GetBinWidth(i))
        h_cross_section.SetBinError(i, h_cross_section.GetBinError(i) / h_cross_section.GetBinWidth(i))

    outfile = ROOT.TFile.Open(config['output_file'], 'recreate')
    h_cross_section.Write()
    outfile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('configFileName', metavar='text', default='config.yaml')
    args = parser.parse_args()
    main(args.configFileName)

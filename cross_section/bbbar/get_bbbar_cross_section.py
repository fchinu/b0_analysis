import argparse
import yaml
import ROOT

def get_bbbar_cross_section(config_file_name):
    '''
    Compute the cross section using the provided configuration file.

    Parameters:
    - config_file_name (str): The path to the configuration file.

    Returns:
    None
    '''

    with open(config_file_name, 'r', encoding='utf-8') as yml_config_file:
        config = yaml.load(yml_config_file, yaml.FullLoader)


    with ROOT.TFile.Open(config['int_cross_sec']) as f:
        h_stat_int_bhadr_cross_sec = f.Get('integrated/h_stat_int')
        h_syst_int_bhadr_cross_sec = f.Get('integrated/h_syst_int_tot')
        g_extrap_int_bhadr_cross_sec = f.Get('integrated/g_extrap')
    
    central_ff = config['FF'][config['FF']['central']]
    rapidity_corr = config['rapidity_correction']['value']
    scale_factor = rapidity_corr/central_ff

    h_stat_bbbar = h_stat_int_bhadr_cross_sec.Clone('h_stat_bbbar')
    h_stat_bbbar.Scale(scale_factor)
    h_syst_bbbar = h_syst_int_bhadr_cross_sec.Clone('h_syst_bbbar')
    h_syst_bbbar.Scale(scale_factor)
    g_extrap_bbbar = g_extrap_int_bhadr_cross_sec.Clone('g_extrap_bbbar')
    g_extrap_bbbar.Scale(scale_factor)
    g_ff_bbbar = g_extrap_bbbar.Clone('g_ff_bbbar')
    g_extrap_ff = g_extrap_bbbar.Clone('g_extrap_ff')

    for i_pt in range(g_ff_bbbar.GetN()):
        differences = [
                h_stat_int_bhadr_cross_sec.GetBinContent(i_pt+1)*rapidity_corr/ff - h_stat_bbbar.GetBinContent(i_pt+1)
                    for ff_name, ff in config['FF'].items() if ff_name != 'central' and ff_name != config['FF']['central']
            ]
        unc_upper = max(
            [diff for diff in differences if diff > 0.], default=0
        )
        unc_lower = max(
            [-diff for diff in differences if diff < 0], default=0
        )
        g_ff_bbbar.SetPointError(i_pt, g_ff_bbbar.GetEXlow()[i_pt], g_ff_bbbar.GetEXhigh()[i_pt], unc_lower, unc_upper)
        extrap_ff_unc_low = ((unc_lower/h_stat_bbbar.GetBinContent(i_pt+1))**2 + (g_extrap_bbbar.GetEYlow()[i_pt]/h_stat_bbbar.GetBinContent(i_pt+1))**2)**0.5 * h_stat_bbbar.GetBinContent(i_pt+1)
        extrap_ff_unc_high = ((unc_upper/h_stat_bbbar.GetBinContent(i_pt+1))**2 + (g_extrap_bbbar.GetEYhigh()[i_pt]/h_stat_bbbar.GetBinContent(i_pt+1))**2)**0.5 * h_stat_bbbar.GetBinContent(i_pt+1)
        g_extrap_ff.SetPointError(i_pt, g_extrap_ff.GetEXlow()[i_pt], g_extrap_ff.GetEXhigh()[i_pt], extrap_ff_unc_low, extrap_ff_unc_high)
    
    with ROOT.TFile.Open(config['output'], 'recreate') as f:
        h_stat_bbbar.Write()
        h_syst_bbbar.Write()
        g_extrap_bbbar.Write()
        g_ff_bbbar.Write()
        g_extrap_ff.Write()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Script to compute the integrated cross section')
    parser.add_argument('config', metavar='text', default='config.yaml')
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)

    get_bbbar_cross_section(args.config)
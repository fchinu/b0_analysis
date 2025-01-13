import argparse
import yaml
import ROOT

def get_cross_sec_with_syst(config_file_name):
    with open(config_file_name, 'r') as f:
        config = yaml.safe_load(f)
    with ROOT.TFile.Open(config['inputs']['cross_section_file']) as f:
        h_cross_section = f.Get('h_cross_section')
        h_lumi_before_bc = f.Get('h_lumi_before_bc')
        h_lumi_after_bc = f.Get('h_lumi_after_bc')

    # Get systematics histograms (pt dependent)
    h_systs_no_br_no_lumi = {}
    h_systs_rel_no_br_no_lumi = {}
    for syst_name, file_name in config['inputs']['syst_files'].items():
        with ROOT.TFile.Open(file_name) as f:
            h_systs_rel_no_br_no_lumi[syst_name] = f.Get('assigned_syst')
            h_systs_no_br_no_lumi[syst_name] = f.Get('assigned_syst')
            h_systs_no_br_no_lumi[syst_name].Multiply(h_cross_section)
    
    # We don't separate tracking syst from the rest
    h_tracking_syst = h_cross_section.Clone()
    h_tracking_syst.Scale(config["tracking"])
    h_tracking_syst_rel = h_cross_section.Clone()
    h_systs_no_br_no_lumi["tracking"] = h_tracking_syst
    for i in range(1, h_cross_section.GetNbinsX()+1):
        h_tracking_syst_rel.SetBinContent(i, config["tracking"])
    h_systs_rel_no_br_no_lumi["tracking"] = h_tracking_syst_rel

    # Get systematics histograms (pt independent)
    h_lumi_syst = h_cross_section.Clone()
    h_lumi_syst.Scale(config["lumi_unc"])
    h_br_syst = h_cross_section.Clone()
    br_unc = 0
    for br in config["br"]:
        br_unc += (br['unc']/br['value'])**2
    br_unc = br_unc**0.5
    h_br_syst.Scale(br_unc)

    h_lumi_syst_rel = h_cross_section.Clone()
    h_br_syst_rel = h_cross_section.Clone()
    for i in range(1, h_cross_section.GetNbinsX()+1):
        h_lumi_syst_rel.SetBinContent(i, config["lumi_unc"])
        h_br_syst_rel.SetBinContent(i, br_unc)
    
    # Initialize histograms for total systematics
    h_total_syst_no_br_no_lumi = h_cross_section.Clone()
    h_total_syst_no_br_no_lumi.Multiply(h_cross_section)
    h_total_syst_w_lumi = h_cross_section.Clone()
    h_total_syst_w_lumi.Multiply(h_cross_section)
    h_total_syst_w_br = h_cross_section.Clone()
    h_total_syst_w_br.Multiply(h_cross_section)
    h_total_syst = h_cross_section.Clone()
    h_total_syst.Multiply(h_cross_section)

    # Initialize histograms for total relative systematics
    h_total_syst_rel_no_br_no_lumi = h_cross_section.Clone()
    h_total_syst_rel_w_lumi = h_cross_section.Clone()
    h_total_syst_rel_w_br = h_cross_section.Clone()
    h_total_syst_rel = h_cross_section.Clone()

    # Initialize histograms for cross section with systematic uncertainties
    h_cross_sec_syst_no_br_no_lumi = h_cross_section.Clone()
    h_cross_sec_syst_no_br_no_lumi.SetName('h_syst_no_br_no_lumi')
    h_cross_sec_syst_w_lumi = h_cross_section.Clone()
    h_cross_sec_syst_w_lumi.SetTitle('h_syst_w_lumi')
    h_cross_sec_syst_w_br = h_cross_section.Clone()
    h_cross_sec_syst_w_br.SetTitle('h_syst_w_br')
    h_cross_sec_syst = h_cross_section.Clone()
    h_cross_sec_syst.SetTitle('h_syst')

    # Evaluate total systematics
    for i in range(1, h_cross_section.GetNbinsX()+1):
        unc = 0
        for syst_name, h in h_systs_no_br_no_lumi.items():
            unc += h.GetBinContent(i)**2
        h_total_syst_no_br_no_lumi.SetBinContent(i, unc**0.5)
        h_total_syst_rel_no_br_no_lumi.SetBinContent(i, unc**0.5/h_cross_section.GetBinContent(i))
        unc_w_lumi = unc + h_lumi_syst.GetBinContent(i)**2
        h_total_syst_w_lumi.SetBinContent(i, unc_w_lumi**0.5)
        h_total_syst_rel_w_lumi.SetBinContent(i, unc_w_lumi**0.5/h_cross_section.GetBinContent(i))
        unc_w_br = unc + h_br_syst.GetBinContent(i)**2
        h_total_syst_w_br.SetBinContent(i, unc_w_br**0.5)
        h_total_syst_rel_w_br.SetBinContent(i, unc_w_br**0.5/h_cross_section.GetBinContent(i))
        unc += h_lumi_syst.GetBinContent(i)**2 + h_br_syst.GetBinContent(i)**2
        h_total_syst.SetBinContent(i, unc**0.5)
        h_total_syst_rel.SetBinContent(i, unc**0.5/h_cross_section.GetBinContent(i))
        

    for i in range(1, h_cross_section.GetNbinsX()+1):
        h_cross_sec_syst_no_br_no_lumi.SetBinError(i, h_total_syst_no_br_no_lumi.GetBinContent(i))
        h_cross_sec_syst_w_lumi.SetBinError(i, h_total_syst_w_lumi.GetBinContent(i))
        h_cross_sec_syst_w_br.SetBinError(i, h_total_syst_w_br.GetBinContent(i))
        h_cross_sec_syst.SetBinError(i, h_total_syst.GetBinContent(i))

        h_total_syst_rel_no_br_no_lumi.SetBinError(i, 0)
        h_total_syst_rel_w_lumi.SetBinError(i, 0)
        h_total_syst_rel_w_br.SetBinError(i, 0)
        h_total_syst_rel.SetBinError(i, 0)

        h_total_syst_no_br_no_lumi.SetBinError(i, 0)
        h_total_syst_w_lumi.SetBinError(i, 0)
        h_total_syst_w_br.SetBinError(i, 0)
        h_total_syst.SetBinError(i, 0)
        for syst_name, h in h_systs_rel_no_br_no_lumi.items():
            h.SetBinError(i, 0)
        for syst_name, h in h_systs_no_br_no_lumi.items():
            h.SetBinError(i, 0)
        h_br_syst.SetBinError(i, 0)
        h_lumi_syst.SetBinError(i, 0)
        h_br_syst_rel.SetBinError(i, 0)
        h_lumi_syst_rel.SetBinError(i, 0)

    with ROOT.TFile.Open(config['output_name'], 'recreate') as f:
        h_cross_section.Write('h_stat')
        h_cross_sec_syst_no_br_no_lumi.Write("h_syst_no_br_no_lumi")
        h_cross_sec_syst_w_lumi.Write('h_syst_w_lumi')
        h_cross_sec_syst_w_br.Write('h_syst_w_br')
        h_cross_sec_syst.Write('h_syst')
        h_total_syst_no_br_no_lumi.Write('h_total_syst_no_br_no_lumi')
        h_total_syst_w_lumi.Write('h_total_syst_w_lumi')
        h_total_syst_w_br.Write('h_total_syst_w_br')
        h_total_syst.Write('h_total_syst')
        h_total_syst_rel_no_br_no_lumi.Write('h_total_syst_rel_no_br_no_lumi')
        h_total_syst_rel_w_lumi.Write('h_total_syst_rel_w_lumi')
        h_total_syst_rel_w_br.Write('h_total_syst_rel_w_br')
        h_total_syst_rel.Write('h_total_syst_rel')
        for syst_name, h in h_systs_no_br_no_lumi.items():
            h.Write(syst_name)
        for syst_name, h in h_systs_rel_no_br_no_lumi.items():
            h.Write(syst_name+'_rel')

        h_lumi_syst.Write('lumi')
        h_br_syst.Write('br')
    
        h_lumi_syst_rel.Write('lumi_rel')
        h_br_syst_rel.Write('br_rel')

        h_lumi_before_bc.Write('h_lumi_before_bc')
        h_lumi_after_bc.Write('h_lumi_after_bc')

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Get cross section with systematics')
    parser.add_argument('config_file', help='Configuration file')
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)

    get_cross_sec_with_syst(args.config_file)

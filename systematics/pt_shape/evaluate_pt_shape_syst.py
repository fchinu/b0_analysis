"""
Script to evaluate pt-shape systematic uncertainty for B0 mesons
"""

import argparse
import pandas as pd
import numpy as np
import yaml
import ROOT

def set_style(histo, color):
    """
    """
    histo.SetDirectory(0)
    histo.SetLineWidth(2)
    histo.SetLineColor(color)
    histo.SetMarkerColor(color)
    histo.SetMarkerStyle(ROOT.kFullCircle)


def evaluate_systematics(infile_fonll, infile_gen, infiles_reco, cutset):
    """
    """

    ROOT.gStyle.SetPadLeftMargin(0.2)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetTitleOffset(2.0, "y")
    ROOT.gStyle.SetTitleOffset(1.4, "x")
    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.045, "xy")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gROOT.SetBatch(True)

    # first we get the pT shapes
    # fonll
    col_names = ["pt", "central", "min", "max", "min_sc", "max_sc",
                 "min_mass", "max_mass", "min_pdf", "max_pdf", "fr_dot5_dot5", "fr_2_2",
                 "fr_2_1", "fr_1_2", "fr_1_dot5", "fr_dot5_1"]
    df = pd.read_csv(infile_fonll, names=col_names, comment="#", sep=" ")

    hist_pt_fonll_cent = ROOT.TH1D("hist_pt_fonll_cent",
                                   ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} normalised",
                                   500, 0., 50.)
    hist_pt_fonll_min = ROOT.TH1D("hist_pt_fonll_min",
                                  ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} normalised",
                                  500, 0., 50.)
    hist_pt_fonll_max = ROOT.TH1D("hist_pt_fonll_max",
                                  ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} normalised",
                                  500, 0., 50.)
    set_style(hist_pt_fonll_cent, ROOT.kAzure+4)
    set_style(hist_pt_fonll_min, ROOT.kAzure+2)
    set_style(hist_pt_fonll_max, ROOT.kBlue+2)

    for ipt, (cent, minval, maxval) in enumerate(
        zip(df["central"].to_numpy(), df["min"].to_numpy(), df["max"].to_numpy())):
        hist_pt_fonll_cent.SetBinContent(ipt+1, cent)
        hist_pt_fonll_min.SetBinContent(ipt+1, minval)
        hist_pt_fonll_max.SetBinContent(ipt+1, maxval)
        hist_pt_fonll_cent.SetBinError(ipt+1, 1.e-20)
        hist_pt_fonll_min.SetBinError(ipt+1, 1.e-20)
        hist_pt_fonll_max.SetBinError(ipt+1, 1.e-20)

    hist_pt_fonll_cent.Scale(1./hist_pt_fonll_cent.Integral())
    hist_pt_fonll_min.Scale(1./hist_pt_fonll_min.Integral())
    hist_pt_fonll_max.Scale(1./hist_pt_fonll_max.Integral())

    # mc gen
    infile_mcgen = ROOT.TFile.Open(infile_gen)
    sparse_gen = infile_mcgen.Get("hf-task-b0-reduced/hPtYGenSig")
    sparse_gen.GetAxis(1).SetRangeUser(-0.499999, 0.499999)
    hist_pt_gen = sparse_gen.Projection(0)
    set_style(hist_pt_gen, ROOT.kRed+1)
    infile_mcgen.Close()
    hist_pt_gen_norm = hist_pt_gen.Clone("hist_pt_gen_norm")
    hist_pt_gen_norm.Scale(1./hist_pt_gen.Integral())

    hist_pt_weights_cent = hist_pt_fonll_cent.Clone("hist_pt_weights_cent")
    hist_pt_weights_cent.Divide(hist_pt_gen_norm)
    hist_pt_weights_min = hist_pt_fonll_min.Clone("hist_pt_weights_min")
    hist_pt_weights_min.Divide(hist_pt_gen_norm)
    hist_pt_weights_max = hist_pt_fonll_max.Clone("hist_pt_weights_max")
    hist_pt_weights_max.Divide(hist_pt_gen_norm)

    # load reco
    df = pd.DataFrame()
    for infile in infiles_reco:
        df = pd.concat([df, pd.read_parquet(infile)])
    # apply selections

    with open(cutset, "r") as yml_cfg:  # pylint: disable=unspecified-encoding
        cfg = yaml.load(yml_cfg, yaml.FullLoader)

    string_selection = ""
    for ipt, (pt_min, pt_max, cut) in enumerate(zip(cfg["pt"]["mins"],
                                                    cfg["pt"]["maxs"],
                                                    cfg["ML_output"]["mins"])): 

        if ipt < len(cfg["pt"]["mins"])-1:
            string_selection += f"({pt_min} < fPt < {pt_max} and ML_output > {cut}) or "
        else:
            string_selection += f"({pt_min} < fPt < {pt_max} and ML_output > {cut})"

    df_sel = df.query(string_selection)

    hist_pt_reco = hist_pt_gen.Clone("hist_reco")
    hist_pt_reco.Reset()
    for pt in df_sel["fPt"].to_numpy():
        hist_pt_reco.Fill(pt)

    hist_pt_reco_fonll_cent = hist_pt_reco.Clone("hist_pt_reco_fonll_cent")
    hist_pt_reco_fonll_min = hist_pt_reco.Clone("hist_pt_reco_fonll_min")
    hist_pt_reco_fonll_max = hist_pt_reco.Clone("hist_pt_reco_fonll_max")
    hist_pt_gen_fonll_cent = hist_pt_gen.Clone("hist_pt_gen_fonll_cent")
    hist_pt_gen_fonll_min = hist_pt_gen.Clone("hist_pt_gen_fonll_min")
    hist_pt_gen_fonll_max = hist_pt_gen.Clone("hist_pt_gen_fonll_max")
    hist_pt_reco_fonll_cent.Sumw2()
    hist_pt_reco_fonll_min.Sumw2()
    hist_pt_reco_fonll_max.Sumw2()
    hist_pt_gen_fonll_cent.Sumw2()
    hist_pt_gen_fonll_min.Sumw2()
    hist_pt_gen_fonll_max.Sumw2()
    set_style(hist_pt_reco_fonll_cent, ROOT.kAzure+4)
    set_style(hist_pt_reco_fonll_min, ROOT.kAzure+2)
    set_style(hist_pt_reco_fonll_max, ROOT.kBlue+2)

    for ipt in range(1, hist_pt_reco.GetNbinsX()+1):
        hist_pt_reco_fonll_cent.SetBinContent(
            ipt, hist_pt_reco_fonll_cent.GetBinContent(ipt) * hist_pt_weights_cent.GetBinContent(ipt))
        hist_pt_reco_fonll_min.SetBinContent(
            ipt, hist_pt_reco_fonll_min.GetBinContent(ipt) * hist_pt_weights_min.GetBinContent(ipt))
        hist_pt_reco_fonll_max.SetBinContent(
            ipt, hist_pt_reco_fonll_max.GetBinContent(ipt) * hist_pt_weights_max.GetBinContent(ipt))
        hist_pt_gen_fonll_cent.SetBinContent(
            ipt, hist_pt_gen_fonll_cent.GetBinContent(ipt) * hist_pt_weights_cent.GetBinContent(ipt))
        hist_pt_gen_fonll_min.SetBinContent(
            ipt, hist_pt_gen_fonll_min.GetBinContent(ipt) * hist_pt_weights_min.GetBinContent(ipt))
        hist_pt_gen_fonll_max.SetBinContent(
            ipt, hist_pt_gen_fonll_max.GetBinContent(ipt) * hist_pt_weights_max.GetBinContent(ipt))

    pt_array = np.array([1, 2, 4, 6, 8, 10, 14, 23.5], dtype=np.float64)
    hist_pt_reco = hist_pt_reco.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_reco_fonll_cent = hist_pt_reco_fonll_cent.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_reco_fonll_min = hist_pt_reco_fonll_min.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_reco_fonll_max = hist_pt_reco_fonll_max.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_gen = hist_pt_gen.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_gen_fonll_cent = hist_pt_gen_fonll_cent.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_gen_fonll_min = hist_pt_gen_fonll_min.Rebin(len(pt_array)-1, "", pt_array)
    hist_pt_gen_fonll_max = hist_pt_gen_fonll_max.Rebin(len(pt_array)-1, "", pt_array)

    hist_eff_noweight = hist_pt_reco.Clone("hist_pt_reco")
    hist_eff_noweight.Divide(hist_pt_reco, hist_pt_gen, 1., 1., "B")
    hist_eff_fonll_cent = hist_pt_reco_fonll_cent.Clone("hist_pt_reco_fonll_cent")
    hist_eff_fonll_cent.Divide(hist_pt_reco_fonll_cent, hist_pt_gen_fonll_cent, 1., 1., "B")
    hist_eff_fonll_min = hist_pt_reco_fonll_min.Clone("hist_pt_reco_fonll_min")
    hist_eff_fonll_min.Divide(hist_pt_reco_fonll_min, hist_pt_gen_fonll_min, 1., 1., "B")
    hist_eff_fonll_max = hist_pt_reco_fonll_max.Clone("hist_pt_reco_fonll_max")
    hist_eff_fonll_max.Divide(hist_pt_reco_fonll_max, hist_pt_gen_fonll_max, 1., 1., "B")

    hist_effratio_fonll_cent = hist_eff_fonll_cent.Clone("hist_effratio_fonll_cent")
    hist_effratio_fonll_cent.Divide(hist_eff_noweight)
    hist_effratio_fonll_min = hist_eff_fonll_min.Clone("hist_effratio_fonll_min")
    hist_effratio_fonll_min.Divide(hist_eff_noweight)
    hist_effratio_fonll_max = hist_eff_fonll_max.Clone("hist_effratio_fonll_max")
    hist_effratio_fonll_max.Divide(hist_eff_noweight)

    for ipt in range(1, hist_effratio_fonll_cent.GetNbinsX()+1):
        hist_effratio_fonll_cent.SetBinError(ipt, 1.e-20)
        hist_effratio_fonll_min.SetBinError(ipt, 1.e-20)
        hist_effratio_fonll_max.SetBinError(ipt, 1.e-20)

    leg = ROOT.TLegend(0.25, 0.2, 0.55, 0.5)
    leg.SetTextSize(0.045)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hist_pt_gen_norm, "MC gen", "pl")
    leg.AddEntry(hist_pt_fonll_cent, "FONLL (central)", "pl")
    leg.AddEntry(hist_pt_fonll_min, "FONLL (min)", "pl")
    leg.AddEntry(hist_pt_fonll_max, "FONLL (max)", "pl")

    canv = ROOT.TCanvas("canv", "", 1000, 1000)
    canv.Divide(2, 2)
    frame = canv.cd(1).DrawFrame(
        0.1, 1.e-5, 50., 1.e-1,
        ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T}d#it{y} normalised"
    )
    canv.cd(1).SetLogx()
    canv.cd(1).SetLogy()
    frame.GetXaxis().SetMoreLogLabels()
    frame.GetYaxis().SetMoreLogLabels()
    hist_pt_gen_norm.DrawCopy("same")
    hist_pt_fonll_cent.DrawCopy("same")
    hist_pt_fonll_min.DrawCopy("same")
    hist_pt_fonll_max.DrawCopy("same")
    leg.Draw()
    frame = canv.cd(2).DrawFrame(
        0.1, 1., 50., 1.,
        ";#it{p}_{T} (GeV/#it{c});FONLL / MC gen"
    )
    canv.cd(2).SetLogx()
    frame.GetXaxis().SetMoreLogLabels()
    frame.GetYaxis().SetDecimals()
    hist_pt_weights_cent.DrawCopy("same")
    hist_pt_weights_min.DrawCopy("same")
    hist_pt_weights_max.DrawCopy("same")
    frame = canv.cd(3).DrawFrame(
        0., 1.e-3, 25., 1.,
        ";#it{p}_{T} (GeV/#it{c});efficiency #times acceptance"
    )
    frame.GetYaxis().SetMoreLogLabels()
    canv.cd(3).SetLogy()
    hist_eff_noweight.DrawCopy("same")
    hist_eff_fonll_cent.DrawCopy("same")
    hist_eff_fonll_min.DrawCopy("same")
    hist_eff_fonll_max.DrawCopy("same")

    frame = canv.cd(4).DrawFrame(
        0., 0.95, 25., 1.05,
        ";#it{p}_{T} (GeV/#it{c});efficiency #times acceptance ratio"
    )
    hist_effratio_fonll_cent.DrawCopy("esame")
    hist_effratio_fonll_min.DrawCopy("esame")
    hist_effratio_fonll_max.DrawCopy("esame")

    canv.SaveAs("pt_shape_syst.pdf")

    outfile = ROOT.TFile("pt_shape_syst.root", "recreate")
    hist_pt_gen_norm.Write()
    hist_pt_fonll_cent.Write()
    hist_pt_fonll_min.Write()
    hist_pt_fonll_max.Write()
    hist_pt_weights_cent.Write()
    hist_pt_weights_min.Write()
    hist_pt_weights_max.Write()
    hist_eff_noweight.Write()
    hist_eff_fonll_cent.Write()
    hist_eff_fonll_min.Write()
    hist_eff_fonll_max.Write()
    hist_effratio_fonll_cent.Write()
    hist_effratio_fonll_min.Write()
    hist_effratio_fonll_max.Write()
    canv.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile_fonll", "-f", metavar="text",
                        default="fonll_shape.txt",
                        help="fonll file", required=False)
    parser.add_argument("--infile_gen", "-g", metavar="text",
                        default="/home/fgrosa/b0_analysis/derived_data_analysis/AnalysisResultsTask_LHC24i4_1.root",
                        help="MC gen file", required=False)
    parser.add_argument("--infile_reco", "-r", nargs="+",
                        default=["LHC24i4_B0ToDPi_pT_1_2_ModelApplied.parquet.gzip",
                                 "LHC24i4_B0ToDPi_pT_2_6_ModelApplied.parquet.gzip",
                                 "LHC24i4_B0ToDPi_pT_6_14_ModelApplied.parquet.gzip",
                                 "LHC24i4_B0ToDPi_pT_14_100_ModelApplied.parquet.gzip"],
                        help="MC reco files with model applied", required=False)
    parser.add_argument("--cutset", "-c", metavar="text", default="cutset.yml",
                        help="config file with cuts", required=False)
    args = parser.parse_args()

    evaluate_systematics(args.infile_fonll, args.infile_gen, args.infile_reco, args.cutset)

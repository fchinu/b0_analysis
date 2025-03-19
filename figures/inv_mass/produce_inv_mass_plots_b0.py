"""

"""
import argparse
import ROOT
import sys
sys.path.append("../../utils")
from style_formatter import root_colors_from_matplotlib_colormap
ci = ROOT.TColor.GetFreeColorIndex()

def set_data_style(hist, isnorm=False):
    """

    """

    hist.SetDirectory(0)
    hist.Sumw2()
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerSize(0.8)
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlack)
    if isnorm:
        hist.GetYaxis().SetTitle(
            f"Normalised counts per {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
        hist.GetYaxis().SetRangeUser(0., hist.GetMaximum() * 1.2)
        hist.SetMarkerSize(1.6)
    else:
        hist.GetYaxis().SetTitle(f"Counts per {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
        hist.GetYaxis().SetRangeUser(0.1, hist.GetMaximum() * 1.2)
    hist.GetXaxis().SetTitle("#it{M}(D^{#mp}#pi^{#pm}) (GeV/#it{c}^{2})")

    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetNdivisions(508)
    hist.GetYaxis().SetNdivisions(508)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetDecimals()
    hist.GetYaxis().SetMaxDigits(3)
    hist.GetXaxis().SetRangeUser(4.9, 5.66)


def get_function_fromhist(hist, func_type="totfunc", norm_fact = 1., colors=None):
    """

    """
    hist.Scale(norm_fact)
    hist.SetFillStyle(1001)
    hist.SetLineWidth(0)

    spline = ROOT.TSpline3(hist)
    spline.SetLineWidth(3)
    if func_type == "totfunc":
        spline.SetLineColor(ROOT.kBlue + 2)
        hist.SetFillColorAlpha(ROOT.kBlue + 2, 0.5)
    elif func_type == "bkg":
        spline.SetLineColor(ROOT.kRed + 1)
        hist.SetLineColor(0)
        hist.SetFillColor(10)
        hist.SetLineStyle(9)
        spline.SetLineStyle(9)
    elif "bkg_corr" in func_type:
        colors_map = {
            0: 3,
            1: 6,
            2: 7,
            3: 1,
            4: 4,
            5: 2
        }
        i_bkg = int(func_type.split("_")[-1])
        spline.SetLineColor(colors[colors_map[i_bkg]*2 + 1])
        hist.SetFillColor(colors[colors_map[i_bkg]*2 + 1])
        spline.SetFillColor(colors[colors_map[i_bkg]*2 + 1])
        spline.SetFillStyle(1000)
        spline.SetLineStyle(i_bkg+2)
    elif func_type == "signal":
        spline.SetLineColor(ROOT.kAzure + 4)
        hist.SetFillColorAlpha(ROOT.kAzure + 4, 0.5)

    return spline

def plot(infile_name, colors, version, pt_min=None, pt_max=None):
    """

    """
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadTopMargin(0.065)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.TGaxis.SetMaxDigits(2)

    if pt_min is None or pt_max is None:
        pt_suffix = "ptint"
    else:
        pt_suffix = f"pt{pt_min}_{pt_max}"

    infile = ROOT.TFile.Open(infile_name)
    hist_mass = infile.Get(f"hdata_{pt_suffix}")
    set_data_style(hist_mass)

    norm = 1./hist_mass.Integral()
    hist_mass_norm = hist_mass.Clone()
    hist_mass_norm.Scale(norm)
    set_data_style(hist_mass_norm, True)
    hist_signal = infile.Get(f"signal_0_{pt_suffix}")
    hist_signal_norm = hist_signal.Clone("hist_signal")
    hist_bkg = infile.Get(f"bkg_6_{pt_suffix}")
    hist_bkg_corr_0 = infile.Get(f"bkg_0_{pt_suffix}")
    hist_bkg_corr_1 = infile.Get(f"bkg_1_{pt_suffix}")
    hist_bkg_corr_2 = infile.Get(f"bkg_2_{pt_suffix}")
    hist_bkg_corr_3 = infile.Get(f"bkg_3_{pt_suffix}")
    hist_bkg_corr_4 = infile.Get(f"bkg_4_{pt_suffix}")
    hist_bkg_corr_5 = infile.Get(f"bkg_5_{pt_suffix}")
    func_bkg = get_function_fromhist(hist_bkg, "bkg")
    func_bkg_corr_0 = get_function_fromhist(hist_bkg_corr_0, "bkg_corr_0", 1., colors)
    func_bkg_corr_1 = get_function_fromhist(hist_bkg_corr_1, "bkg_corr_1", 1., colors)
    func_bkg_corr_2 = get_function_fromhist(hist_bkg_corr_2, "bkg_corr_2", 1., colors)
    func_bkg_corr_3 = get_function_fromhist(hist_bkg_corr_3, "bkg_corr_3", 1., colors)
    func_bkg_corr_4 = get_function_fromhist(hist_bkg_corr_4, "bkg_corr_4", 1., colors)
    func_bkg_corr_5 = get_function_fromhist(hist_bkg_corr_5, "bkg_corr_5", 1., colors)
    if version == "2":
        corr_bkg_stack = ROOT.THStack("corr_bkg_stack", "")
        corr_bkg_stack.Add(hist_bkg_corr_1)
        corr_bkg_stack.Add(hist_bkg_corr_2)
        corr_bkg_stack.Add(hist_bkg_corr_3)
        corr_bkg_stack.Add(hist_bkg_corr_0)
        corr_bkg_stack.Add(hist_bkg_corr_4)
        corr_bkg_stack.Add(hist_bkg_corr_5)
    if version == "3":
        corr_bkg_stack = ROOT.THStack("corr_bkg_stack", "")
        corr_bkg_stack.Add(hist_bkg)
        corr_bkg_stack.Add(hist_bkg_corr_1)
        corr_bkg_stack.Add(hist_bkg_corr_2)
        corr_bkg_stack.Add(hist_bkg_corr_3)
        corr_bkg_stack.Add(hist_bkg_corr_0)
        corr_bkg_stack.Add(hist_bkg_corr_4)
        corr_bkg_stack.Add(hist_bkg_corr_5)

    func_signal = get_function_fromhist(hist_signal, "signal")
    func_totfunc = get_function_fromhist(infile.Get(f"total_func_{pt_suffix}"), "totfunc")

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextColor(ROOT.kBlack)
    lat.SetTextSize(0.045)

    lat_small = ROOT.TLatex()
    lat_small.SetNDC()
    lat_small.SetTextFont(42)
    lat_small.SetTextColor(ROOT.kBlack)
    lat_small.SetTextSize(0.04)

    leg = ROOT.TLegend(0.62, 0.5, 0.92, 0.7)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hist_mass, "Data", "p")
    leg.AddEntry(hist_signal, "#splitline{#lower[0.05]{B^{0}#rightarrow D^{#minus}#pi^{#plus}}}{#lower[-0.05]{and charge conj.}}", "f")
    leg.AddEntry(func_bkg, "Comb. background", "l")
    leg.AddEntry(func_totfunc, "Total fit function", "l")

    leg2 = ROOT.TLegend(0.62, 0.5, 0.92, 0.7)
    leg2.SetTextSize(0.03)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.AddEntry(hist_mass, "Data", "p")
    leg2.AddEntry(func_signal, "#splitline{#lower[0.05]{B^{0}#rightarrow D^{#minus}#pi^{#plus}#rightarrow#pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}}}{#lower[-0.05]{and charge conj.}}", "f")
    leg2.AddEntry(func_bkg, "Comb. background", "l")
    leg2.AddEntry(func_totfunc, "Total fit function", "l")

    if pt_min is None and pt_max is None:
        leg_corr = ROOT.TLegend(0.14, 0.3, 0.35, 0.55)
    elif pt_min==2 and pt_max==4:
        leg_corr = ROOT.TLegend(0.14, 0.2, 0.35, 0.45)
    elif pt_min==10 and pt_max==14:
        leg_corr = ROOT.TLegend(0.14, 0.23, 0.35, 0.48)
    elif pt_min==14 and pt_max==24:
        leg_corr = ROOT.TLegend(0.14, 0.28, 0.35, 0.53)
    leg_corr.SetTextSize(0.03)
    leg_corr.SetBorderSize(0)
    leg_corr.SetFillStyle(0)
    if version == 1:
        leg_corr.AddEntry(func_bkg_corr_5, "B^{0}#rightarrow D^{#minus}K^{+}", "l")
        leg_corr.AddEntry(func_bkg_corr_1, "B^{0}#rightarrow D^{*#minus}#pi^{+}#rightarrow D^{0}#pi^{#minus}#pi^{+}", "l")
        leg_corr.AddEntry(func_bkg_corr_2, "B_{s}^{0}#rightarrow D_{s}^{#minus}#pi^{+}#rightarrow K^{+}K^{#minus}#pi^{+}#pi^{+}", "l")
        leg_corr.AddEntry(func_bkg_corr_3, "#Lambda_{b}^{0}#rightarrow #kern[-0.5]{#Lambda_{c}^{+}}#pi^{#minus}#rightarrow pK^{#minus}#pi^{+}#pi^{#minus}", "l")
        leg_corr.AddEntry(func_bkg_corr_0, "B^{0}#rightarrow D^{*#minus}#pi^{+}#rightarrow D^{#minus}#pi^{+}{#pi^{0},#gamma}", "l")
        leg_corr.AddEntry(func_bkg_corr_4, "B^{0}#rightarrow D^{#minus}#rho^{+}#rightarrow D^{#minus}#pi^{+}{#pi^{0},#gamma}", "l")
    else:
        leg_corr.AddEntry(hist_bkg_corr_5, "B^{0}#rightarrow D^{#minus}K^{+}", "f")
        leg_corr.AddEntry(hist_bkg_corr_1, "B^{0}#rightarrow D^{*#minus}#pi^{+}#rightarrow D^{0}#pi^{#minus}#pi^{+}", "f")
        leg_corr.AddEntry(hist_bkg_corr_2, "B_{s}^{0}#rightarrow D_{s}^{#minus}#pi^{+}#rightarrow K^{+}K^{#minus}#pi^{+}#pi^{+}", "f")
        leg_corr.AddEntry(hist_bkg_corr_3, "#Lambda_{b}^{0}#rightarrow #kern[-0.5]{#Lambda_{c}^{+}}#pi^{#minus}#rightarrow pK^{#minus}#pi^{+}#pi^{#minus}", "f")
        leg_corr.AddEntry(hist_bkg_corr_0, "B^{0}#rightarrow D^{*#minus}#pi^{+}#rightarrow D^{#minus}#pi^{+}{#pi^{0},#gamma}", "f")
        leg_corr.AddEntry(hist_bkg_corr_4, "B^{0}#rightarrow D^{#minus}#rho^{+}#rightarrow D^{#minus}#pi^{+}{#pi^{0},#gamma}", "f")

    canv_masses = ROOT.TCanvas("canv_masses", "", 500, 500)
    hist_mass.DrawCopy()
    if version == "1":
        func_bkg.Draw("lsame")
        func_bkg_corr_0.Draw("lsame")
        func_bkg_corr_1.Draw("lsame")
        func_bkg_corr_2.Draw("lsame")
        func_bkg_corr_3.Draw("lsame")
        func_bkg_corr_4.Draw("lsame")
        func_bkg_corr_5.Draw("lsame")
        func_totfunc.Draw("lsame")
        hist_signal.DrawCopy("histsame")
        func_signal.Draw("lsame")
    elif version == "2":
        func_bkg.Draw("lsame")
        corr_bkg_stack.Draw("NOCLEAR hist same")
        func_totfunc.Draw("lsame")
        hist_signal.DrawCopy("histsame")
        func_signal.Draw("lsame")
    elif version == "3":
        corr_bkg_stack.Draw("NOCLEAR hist same")
        hist_bkg.Draw("histsame")
        func_bkg.Draw("lsame")
        hist_mass.DrawCopy("same")
        func_totfunc.Draw("lsame")
        hist_signal.DrawCopy("histsame")
        func_signal.Draw("lsame")
    lat.DrawLatex(0.595, 0.86, "ALICE Preliminary")
    lat.DrawLatex(0.35, 0.8,
                  "pp,#sqrt{#it{s}} = 13.6 TeV, #font[132]{#it{L}}_{int} = 43 pb^{#minus1}")

    if pt_min is None and pt_max is None:
        lat.DrawLatex(0.56, 0.74, "1 < #kern[-0.3]{#it{p}_{T}} < 23.5 GeV/#it{c}")
    elif pt_min==2 and pt_max==4:
        lat.DrawLatex(0.60, 0.74, f"{pt_min} < #it{{p}}_{{T}} < {pt_max} GeV/#it{{c}}")
    elif pt_min==10 and pt_max==14:
        lat.DrawLatex(0.56, 0.74, f"{pt_min} < #it{{p}}_{{T}} < {pt_max} GeV/#it{{c}}")
    leg.Draw()
    leg2.Draw()
    leg_corr.Draw()

    canv_masses.RedrawAxis()

    canv_masses.SaveAs(f"B0_mass_{pt_suffix}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        default="../b0_analysis/fit/B0_mass.root",
                        help="input root file", required=False)
    parser.add_argument("--palette", "-p", metavar="text",
                        default="tab20",
                        help="matplotlib palette name for corr. bkg", required=False)
    parser.add_argument("--version", "-v", metavar="text",
                        default="1",
                        help="style version", required=False)
    parser.add_argument("--ptmin", type=int,
                        default=None,
                        help="minimum pt", required=False)
    parser.add_argument("--ptmax", type=int,
                        default=None,
                        help="maximum pt", required=False)
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)
    COLORS, _ = root_colors_from_matplotlib_colormap(args.palette)

    plot(args.input_file, COLORS, args.version, args.ptmin, args.ptmax)

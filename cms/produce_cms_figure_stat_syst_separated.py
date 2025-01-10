import pandas as pd
import numpy as np
import ROOT

data = pd.read_csv("cms/cms_bplus_13_tev_table_2.txt", sep=' ', comment='#')
data["syst_tot"] = np.sqrt(data["syst"]**2 + data["lumi_unc"]**2)

pt_edges = np.array(data["pt_min"]).tolist() + [data["pt_max"].iloc[-1]]

g_stat = ROOT.TGraphAsymmErrors(len(pt_edges)-1)
g_syst = ROOT.TGraphAsymmErrors(len(pt_edges)-1)

for i in range(len(pt_edges)-1):
    pt = (data["pt_min"].loc[i] + data["pt_max"].loc[i])/2
    g_stat.SetPoint(i, pt, data["sigma_central"].iloc[i])
    g_stat.SetPointError(
        i, data["pt_max"].loc[i]-pt, pt-data["pt_min"].loc[i],
        data["stat_upper"].iloc[i], data["stat_lower"].iloc[i]
    )
    g_syst.SetPoint(i, pt, data["sigma_central"].iloc[i])
    g_syst.SetPointError(
        i, data["pt_max"].loc[i]-pt, pt-data["pt_min"].loc[i],
        data["syst_tot"].iloc[i], data["syst_tot"].iloc[i]
    )

with ROOT.TFile.Open("cms/cms_bplus_13_tev_table_2.root", "recreate") as f:
    g_stat.Write("g_stat")
    g_syst.Write("g_syst")

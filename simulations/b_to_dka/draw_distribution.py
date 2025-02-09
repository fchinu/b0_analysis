import uproot
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import ROOT

file = ROOT.TFile.Open('simulations/b_to_dka/MDPi_511_from_DKa_decay.root')
tree = file.Get('treeBtoDKa')

# Draw the distributions

c = ROOT.TCanvas('c', 'c', 700, 600)
h_frame = c.DrawFrame(4.9, 0, 5.66, 0.4, ";m_{B} (GeV/c^{2});Counts (a.u.)")

tree.Draw('mB >> h_mb(76,4.9,5.66)', 'pTKa<3', 'goff')
h_mb = ROOT.gDirectory.Get('h_mb')

tree.Draw('mBcorr >> h_mbcorr(76,4.9,5.66)', 'pTKa<3', 'goff')
h_mbcorr = ROOT.gDirectory.Get('h_mbcorr')

tree.Draw('mBSmeared >> h_mb_smeared(76,4.9,5.66)', 'pTKa<3', 'goff')
h_mb_smeared = ROOT.gDirectory.Get('h_mb_smeared')

tree.Draw('mBcorrSmeared >> h_mbcorr_smeared(76,4.9,5.66)', 'pTKa<3', 'goff')
h_mbcorr_smeared = ROOT.gDirectory.Get('h_mbcorr_smeared')

h_mb.SetLineColor(ROOT.kRed)
h_mb.SetLineWidth(2)
h_mb.SetMarkerStyle(ROOT.kFullCircle)
h_mb.SetMarkerSize(1.5)
h_mb.SetMarkerColor(ROOT.kRed)
h_mbcorr.SetLineColor(ROOT.kAzure-3)
h_mbcorr.SetLineWidth(2)
h_mbcorr.SetMarkerStyle(ROOT.kFullDiamond)
h_mbcorr.SetMarkerSize(2)
h_mbcorr.SetMarkerColor(ROOT.kAzure-3)

h_mb_smeared.SetLineColor(ROOT.kRed)
h_mb_smeared.SetLineWidth(0)
h_mb_smeared.SetMarkerStyle(ROOT.kFullCircle)
h_mb_smeared.SetMarkerSize(1.5)
h_mb_smeared.SetMarkerColor(ROOT.kRed)
h_mb_smeared.SetFillColorAlpha(ROOT.kRed, 0.5)
h_mb_smeared.SetFillStyle(1001)
h_mbcorr_smeared.SetLineColor(ROOT.kAzure-3)
h_mbcorr_smeared.SetLineWidth(0)
h_mbcorr_smeared.SetMarkerStyle(ROOT.kFullDiamond)
h_mbcorr_smeared.SetMarkerSize(2)
h_mbcorr_smeared.SetMarkerColor(ROOT.kAzure-3)
h_mbcorr_smeared.SetFillColorAlpha(ROOT.kAzure-3, 0.5)
h_mbcorr_smeared.SetFillStyle(1001)

h_mb.DrawNormalized('same')
h_mbcorr.DrawNormalized('same')
h_mb_smeared.DrawNormalized('same')
h_mbcorr_smeared.DrawNormalized('same')

leg = ROOT.TLegend(0.6, 0.6, 0.8, 0.8)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.AddEntry(h_mb, 'm_{B} (K as #pi)', 'lp')
leg.AddEntry(h_mb_smeared, 'm_{B} (K as #pi), smeared', 'f')
leg.AddEntry(h_mbcorr, 'm_{B}', 'lp')
leg.AddEntry(h_mbcorr_smeared, 'm_{B}, smeared', 'f')
leg.Draw()

c.SaveAs('simulations/b_to_dka/mass_B.png')
file.Close()

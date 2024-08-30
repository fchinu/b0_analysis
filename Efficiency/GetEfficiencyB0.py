'''
Script for the computation of the B0 meson efficiency
run: python GetEfficiencyB0.py configFileName.yml cutSetFileName.yml
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
from AnalysisUtils import EvaluateEfficiencyFromHistos

def ComputeEfficiency(configFileName, cutSetFileName):
    """
    Compute the efficiency of a B0 particle based on the given configuration and cut set.

    Args:
        configFileName (str): The file name of the configuration file.
        cutSetFileName (str): The file name of the cut set file.

    Returns:
        None
    """
    
    with open(args.cutSetFileName, 'r') as ymlcutsetFile:
        cutset = yaml.safe_load(ymlcutsetFile)

    with open(args.configFileName, 'r') as ymlconfigFile:
        config = yaml.safe_load(ymlconfigFile)

    # Set style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetTitleOffset(1.2, 'X')
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ptMins = cutset['pt']['mins']
    ptMaxs = cutset['pt']['maxs']
    ptLims = list(ptMins)
    nPtBins = len(ptMins)
    ptLims.append(ptMaxs[-1])

    hRecoIntegrated = ROOT.TH1F('hRecoIntegrated', ';#it{p}_{T} (GeV/#it{c});Reconstructed', nPtBins, np.asarray(ptLims, 'd'))
    hGenIntegrated = ROOT.TH1F('hGenIntegrated', ';#it{p}_{T} (GeV/#it{c});Generated', nPtBins, np.asarray(ptLims, 'd'))

    for iFile, inFileName in enumerate(config['gen']['fileNames']):
        inFile = ROOT.TFile.Open(inFileName)
        hSparseGen = inFile.Get(config['gen']['sparseName'])
        if iFile == 0:
            hGen = hSparseGen.Projection(config['gen']['ptAxis'])
            hGen.SetDirectory(0)
        else:
            hGen.Add(hSparseGen.Projection(config['gen']['ptAxis']))
        inFile.Close()

    hGen.SetName('hGen')
    hReco = hGen.Clone('hReco')
    hReco.Reset()


    for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):

        dfReco = pd.concat([read_parquet_in_batches(parquet, f"{ptMin} < fPt < {ptMax}") for parquet in config['recoFileNames']])
        
        selToApply = ''
        for cutVar in cutset:
            if cutVar == 'pt' or cutVar == 'M':
                continue
            selToApply += f" {cutset[cutVar]['mins'][iPt]} < {cutVar} < {cutset[cutVar]['maxs'][iPt]} and"
        selToApply = selToApply[:-3] # Remove the last 'and'

        dfReco = dfReco.query(selToApply)
        for pt in dfReco['fPt']:
            hReco.Fill(pt)

        # Get the number of generated and reconstructed particles in the given pt range
        nRecoUnc, nGenUnc = (ctypes.c_double() for _ in range(2))
        nReco = hReco.IntegralAndError(hReco.FindBin(ptMin), hReco.FindBin(ptMax)-1, nRecoUnc)
        nGen = hGen.IntegralAndError(hGen.FindBin(ptMin), hGen.FindBin(ptMax)-1, nGenUnc)

        hRecoIntegrated.SetBinContent(iPt+1, nReco)
        hRecoIntegrated.SetBinError(iPt+1, nRecoUnc.value)
        hGenIntegrated.SetBinContent(iPt+1, nGen)
        hGenIntegrated.SetBinError(iPt+1, nGenUnc.value)

    # Compute the efficiency in the given pt ranges
    hEff = EvaluateEfficiencyFromHistos(hGenIntegrated, hRecoIntegrated)
    hEff.SetMarkerStyle(ROOT.kFullCircle)
    hEff.SetMarkerColor(ROOT.kAzure+3)
    hEff.SetMarkerSize(1.5)
    hEff.SetLineColor(ROOT.kAzure+3)
    hEff.SetLineWidth(2)

    # Compute the efficiency in pt bins used for the generated particles
    hEffFineBins = EvaluateEfficiencyFromHistos(hGen, hReco)

    leg = ROOT.TLegend(0.6, 0.2, 0.8, 0.4)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hEff, " B^{0}", "p")

    cEff = ROOT.TCanvas('cEff', '', 800, 800)
    cEff.DrawFrame(0, 1.e-4, ptMaxs[nPtBins-1], 1.,
                ';#it{p}_{T} (GeV/#it{c});Efficiency;')
    cEff.SetLogy()
    hEff.Draw('same')
    leg.Draw()

    outFile = ROOT.TFile(config['outputFileName'], 'recreate')
    hEff.Write()
    hEffFineBins.Write()
    hGen.Write()
    hReco.Write()
    hGenIntegrated.Write()
    hRecoIntegrated.Write()
    outFile.Close()

    outFileNamePDF = config['outputFileName'].replace('.root', '.pdf')
    cEff.SaveAs(outFileNamePDF)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('cutSetFileName', metavar='text', default='config_Cutset_B0.yml')
    parser.add_argument('configFileName', metavar='text', default='config_Efficiency')
    args = parser.parse_args()

    ComputeEfficiency(args.cutSetFileName, args.configFileName)

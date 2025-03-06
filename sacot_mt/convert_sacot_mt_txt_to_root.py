import argparse

import pandas as pd 
import numpy as np
import ROOT 

def convert_to_root(input_file, output_file):
    df = pd.read_csv(input_file, delimiter=r"\s+", comment='#')
    print(df)
    pt_mins = df['pTmin'].values
    pt_maxs = df['pTmax'].values
    pt_bins = list(pt_mins) + [pt_maxs[-1]]
    
    h = ROOT.TH1F('h', 'h', len(pt_bins)-1, np.asarray(pt_bins, "d"))
    for i in range(1, len(pt_bins)):
        h.SetBinContent(i, df['central'].values[i-1])
        h.SetBinError(i, 0)
    
    g = ROOT.TGraphAsymmErrors(h)
    for i in range(1, len(pt_bins)):
        g.SetPointEYlow(i-1, (df['central'].values[i-1] - df['minimum'].values[i-1]))
        g.SetPointEYhigh(i-1, (df['maximum'].values[i-1] - df['central'].values[i-1])) 
   
    h.Scale(0.5)    # 0.5 is the factor to account for the fact that the cross section is for B0 and B0bar
    g.Scale(0.5)    # 0.5 is the factor to account for the fact that the cross section is for B0 and B0bar
    # h.Scale(0.407)  # FF
    # g.Scale(0.407)  # FF

    output_file = ROOT.TFile(output_file, 'RECREATE')
    h.Write("hBhadr")
    g.Write("gBhadr")
    output_file.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert sacot mT cross section to ROOT file')
    parser.add_argument('input_file', type=str, help='Input file')
    parser.add_argument('output_file', type=str, help='Output file')
    args = parser.parse_args()

    convert_to_root(args.input_file, args.output_file)
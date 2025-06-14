import argparse

import pandas as pd 
import numpy as np
import ROOT 

def convert_to_root(input_file, output_file):
    df = pd.read_csv(input_file, delimiter=r"\s+", names=['pT(GeV)', 'min', 'max'])
    df['central'] = (df['min'] + df['max']) / 2
    print(df)
    
    g = ROOT.TGraphAsymmErrors(len(df))
    for i in range(len(df)):
        g.SetPoint(i, df['pT(GeV)'].values[i], df['central'].values[i])
        g.SetPointError(i,
            0.1, 0.1,
            (df['central'].values[i] - df['min'].values[i]) / 2,
            (df['max'].values[i] - df['central'].values[i]) / 2
        )

    output_file = ROOT.TFile(output_file, 'RECREATE')
    g.Write("gBhadr")
    output_file.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert FONLL cross section to ROOT file')
    parser.add_argument('input_file', type=str, help='Input file')
    parser.add_argument('output_file', type=str, help='Output file')
    args = parser.parse_args()

    convert_to_root(args.input_file, args.output_file)
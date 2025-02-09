import argparse

import pandas as pd 
import numpy as np
import ROOT 

def convert_to_root(input_file, output_file):
    df = pd.read_csv(input_file, delimiter=r"\s+", skiprows=14)
    print(df)
    
    g = ROOT.TGraph(len(df))
    for i in range(len(df)):
        g.SetPoint(i, df['pT(GeV)'].values[i], df['B-=B0bar'].values[i])

    output_file = ROOT.TFile(output_file, 'RECREATE')
    g.Write("gBhadr")
    output_file.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert FONLL cross section to ROOT file')
    parser.add_argument('input_file', type=str, help='Input file')
    parser.add_argument('output_file', type=str, help='Output file')
    args = parser.parse_args()

    convert_to_root(args.input_file, args.output_file)
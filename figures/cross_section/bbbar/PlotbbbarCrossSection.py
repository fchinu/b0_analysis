'''
Script for the plots of the ds_bbbar/dy cross section vs other measurements and sqrts
'''

import sys
import ctypes
import numpy as np
import yaml
from ROOT import TFile, TCanvas, TLegend, TLatex, TGraphAsymmErrors, TLine, TColor, TH2F # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kRed, kBlack, kGray, kOrange, kMagenta, kGreen, kFullDiamond, kOpenCircle, kOpenSquare, kOpenDiamond # pylint: disable=import-error,no-name-in-module
from ROOT import kOpenCross, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kFullCircle, kFullCross, kFullDoubleDiamond # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
sys.path.append('../../../utils')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, kRedMy # pylint: disable=import-error,no-name-in-module
from utils.ReadModel import ReadFONLL
from style_formatter import root_colors_from_matplotlib_colormap

colors, col_idx = root_colors_from_matplotlib_colormap("tab10")

# define custom colors
kAzureMy = TColor.GetFreeColorIndex()
cAzureMy = TColor(kAzureMy, 217./255, 229./255, 242./255, 'kAzureMy', 1.0)
# azure a bit darker 
kAzureMyDarker = TColor.GetFreeColorIndex()
cAzureMyDarker = TColor(kAzureMyDarker, 207./255, 219./255, 232./255, 'kAzureMyDarker', 1.0)
kGreenMy = TColor.GetFreeColorIndex()
cGreenMy = TColor(kGreenMy, 179./255, 230./255, 179./255, 'kGreenMy', 1.0)

# useful variables
measToPlot13TeV = ['DzeroLc', 'Lc', 'Dzero', 'Jpsi', 'eePO', 'eePy', 'Dplus', 'Ds']
measToPlotVsSqrt = {'200GeV': ['HFe'],
                    '630GeV': ['B'],
                    '1dot96TeV': ['Jpsi'],
                    '5dot02TeV': ['eePy', 'eePO', 'Jpsi', 'av'],
                    '7TeV': ['eePy', 'eePO', 'Jpsi', 'HFe'],
                    '13TeV': ['eePy', 'eePO', 'Jpsi', 'DzeroLc'],
                    '13dot6TeV': ['B0']}

# compilation plot
yTitles = {'DzeroLc': 'b#rightarrow D + b#rightarrow #Lambda_{c}^{+} |#it{y}|<0.5', 
           'av': 'b#rightarrow D^{0}, #Lambda_{c}^{+} average |#it{y}|<0.5', 'Lc': 'b#rightarrow #Lambda_{c}^{+} |#it{y}|<0.5',
           'Dzero': 'b#rightarrow D^{0} |#it{y}|<0.5', 'Jpsi': 'b#rightarrow J/#psi |#it{y}|<0.9',
           'eePO': 'Dielectron |#it{#eta}_{e}|<0.8 POWHEG', 'eePy': 'Dielectron |#it{#eta}_{e}|<0.8 PYTHA6',
           'Dplus': 'b#rightarrow D^{+} |#it{y}|<0.5', 'Ds': 'b#rightarrow D_{s}^{+} |#it{y}|<0.5',
           'empty': ''}

order = {'DzeroLc': 1, 'Lc': 2, 'Dzero': 3, 'Dplus': 4, 'Ds': 5, 'Jpsi': 6, 'eePO': 7, 'eePy': 8, 'empty': 9, 'av': -1}
# references = {'av': '', 'Ds': '', 'Dplus': '', 'Dzero': '',
#              'Jpsi': 'Preliminary', 'eePO': 'arXiv:2005.11995', 'eePy': 'arXiv:2005.11995', 'empty': ''}
references = {'DzeroLc': '', 'av': '', 'Lc': '', 'Dzero': '',
              'Jpsi': '', 'eePO': '', 'eePy': '', 'empty': '', 'Dplus': '', 'Ds': ''}
refPosX = {'DzeroLc': -20, 'av': -20, 'Lc': -20, 'Dzero': -20,
           'Jpsi': -17, 'eePO': -25.5, 'eePy': -25.5, 'empty': -20,
           'Dplus': -20, 'Ds': -20}

# vs sqrts
energies = {'68dot6GeV': 0.0686, '200GeV': 0.2, '500GeV': 0.5, '630GeV': 0.63, '900GeV': 0.9, '1dot96TeV': 1.96,
            '2dot76TeV': 2.76, '5dot02TeV': 5.02, '5dot36TeV': 5.36, '5dot5TeV': 5.5, '7TeV': 7, '8TeV': 8,
            '8dot16TeV': 8.16, '13TeV': 13, '13dot6TeV': 13.6}

# energy shift for better visibility
shift = {'Jpsi': 1., 'eePO': 1.,
         'eePy': 1., 'HFe': 1., 'B': 1., 'av': 1., 'Dzero': 1., 'Lc': 1., 'DzeroLc': 1., 'B0': 1.}
markerstyles = {'200GeV': {'HFe': kOpenStar}, '630GeV': {'B': kOpenSquare}, '1dot96TeV': {'Jpsi': kOpenDiamond},
                '5dot02TeV': {'eePy': kOpenTriangleUp, 'eePO': kOpenTriangleDown, 'Jpsi': kOpenCircle},
                '7TeV': {'eePy': kOpenTriangleUp, 'eePO': kOpenTriangleDown, 'Jpsi': kOpenCircle, 'HFe': kOpenCross},
                '13TeV': {'eePy': kOpenTriangleUp, 'eePO': kOpenTriangleDown, 'Jpsi': kOpenCircle},
                '13dot6TeV': {'B0': kFullCircle}}

# other measurements
# from HepData when available, hard codded otherwise
# **********************************************************************************************
# PHENIX 200 GeV Phys. Rev. Letters 103 (2009) 082002: HepData (table 2)

# UA1 630 GeV Phys. Lett. B 256 (1991) 121-128 : HepData (table 10)

# CDF 1.96 TeV Phys. Rev. D 71 (2005) 032001: HepData (table 5)

# ALICE 5 TeV (dielectrons PYTHIA) arXiv:2005.11995: 34±4(stat)±2(sys)±2(BR)
# ALICE 5 TeV (dielectrons POWHEG) arXiv:2005.11995: 28±5(stat)±1(sys)±2(BR)
# ALICE 5 TeV (np J/psi): 34.7±5.9(stat)±5.2(syst)±3.0(BR)±0.7(lumi)+1.1-2.3(extr.) mub

# ALICE 7 TeV (HFe) Phys. Lett. B 721 (2013) 13-23: 42.3±3.5(stat)+12.3−11.9(sys)+1.1−1.7(extr)
# ALICE 7 TeV (np J/psi) JHEP 1211 (2012) 065, 2012: HepData (table 6)
# ALICE 7 TeV (dielectrons PYTHIA) JHEP 1809 (2018) 064: 48 ± 8 (stat) ± 12 (syst) ± 3 (BR) = 48 ± 17.7 mub
# ALICE 7 TeV (dielectrons POWHEG) JHEP 1809 (2018) 064: 23 ± 4 (stat) ± 11 (syst) ± 1 (BR) = 23 ± 11.75 mub

# ALICE 13 TeV (dielectrons PYTHIA) Phys. Lett. B 788 (2019) 505: 79±14(stat)±11(sys)
# ALICE 13 TeV (dielectrons POWHEG) Phys. Lett. B 788 (2019) 505: 48±14(stat.)±7(sys)
# ALICE 13 TeV (np J/psi) JHEP 03 (2022) 190: 73.3±6.1(stat)±9.3(sys)+0.8−2.3(extr)
#
# ALICE 13.6 TeV (B0) Preliminary: 82.2±4.6(stat)±12.9(sys)+0.01-12.9(extr)
# **********************************************************************************************

sigma = {'5dot02TeV': {'eePO': 28, 'eePy': 34, 'Jpsi': 34.7},
         '7TeV': {'eePO': 23, 'eePy': 48, 'HFe': 42.3},
         '13TeV': {'eePO': 48, 'eePy': 79, 'Jpsi': 73.3},
         '13dot6TeV': {'B0': 80.6}}

stat = {'5dot02TeV': {'eePO': 5, 'eePy': 4, 'Jpsi': 5.9},
        '7TeV': {'eePO': 4, 'eePy': 8, 'HFe': 3.5},
        '13TeV': {'eePO': 14, 'eePy': 14, 'Jpsi': 6.1},
        '13dot6TeV': {'B0': 4.4}}

systLow = {'5dot02TeV': {'eePO': np.sqrt(1**2+2**2), 'eePy': np.sqrt(2**2+2**2),
                         'Jpsi': np.sqrt(5.2**2+3.0**2+0.7**2)},
           '7TeV': {'eePO': np.sqrt(11**2+1**2), 'eePy': np.sqrt(12**2+3**2), 'HFe': 11.9},
           '13TeV': {'eePO': 7, 'eePy': 11, 'Jpsi': 9.3},
           '13dot6TeV': {'B0': 12.7}}

systHigh = {'5dot02TeV': {'eePO': np.sqrt(1**2+2**2), 'eePy': np.sqrt(2**2+2**2),
                          'Jpsi': np.sqrt(5.2**2+3.0**2+0.7**2)},
            '7TeV': {'eePO': np.sqrt(11**2+1**2), 'eePy': np.sqrt(12**2+3**2), 'HFe': 12.3},
            '13TeV': {'eePO': 7, 'eePy': 11, 'Jpsi': 9.3},
            '13dot6TeV': {'B0': 12.7}}

extrapSystLow = {'5dot02TeV': {'Jpsi': 2.3, 'eePO': 0, 'eePy': 0},
                 '7TeV': {'eePO': 0, 'eePy': 0, 'HFe': 1.7},
                 '13TeV': {'eePO': 0, 'eePy': 0, 'Jpsi': 2.3},
                 '13dot6TeV': {'B0': 12.7}}

extrapSystHigh = {'5dot02TeV': {'Jpsi': 1.1, 'eePO': 0, 'eePy': 0},
                  '7TeV': {'eePO': 0, 'eePy': 0, 'HFe': 1.1},
                  '13TeV': {'eePO': 0, 'eePy': 0, 'Jpsi': 0.8},
                  '13dot6TeV': {'B0': 0.01}}

E, dsdy = ctypes.c_double(), ctypes.c_double()

inFilePHENIX = TFile.Open('ptintcross/HEPData_bbbar_PHENIX-pp200GeV_PRL_103-082002_2009.root')
gTmp = inFilePHENIX.Get('Table 2/Graph1D_y1')
gTmp.GetPoint(0, E, dsdy)
sigma['200GeV'], stat['200GeV'], systLow['200GeV'], systHigh['200GeV'], \
    extrapSystLow['200GeV'], extrapSystHigh['200GeV'] = ({} for _ in range(6))
sigma['200GeV']['HFe'] = dsdy.value
stat['200GeV']['HFe'] = 0.  # info not available
systLow['200GeV']['HFe'] = gTmp.GetErrorYlow(0)
systHigh['200GeV']['HFe'] = gTmp.GetErrorYhigh(0)
extrapSystLow['200GeV']['HFe'] = 0.
extrapSystHigh['200GeV']['HFe'] = 0.
inFilePHENIX.Close()

inFileUA1 = TFile.Open('ptintcross/HEPData_bbbar_UA1_ppbar630GeV_PLB_256_121-128_1991.root')
gTmp = inFileUA1.Get('Table 10/Graph1D_y1')
gTmp.GetPoint(0, E, dsdy)
sigma['630GeV'], stat['630GeV'], systLow['630GeV'], systHigh['630GeV'], \
    extrapSystLow['630GeV'], extrapSystHigh['630GeV'] = ({} for _ in range(6))
sigma['630GeV']['B'] = dsdy.value / 3  # measurement integrated in |y|<1.5
stat['630GeV']['B'] = 0.  # info not available
systLow['630GeV']['B'] = gTmp.GetErrorYlow(
    0) / 3  # measurement integrated in |y|<0.6
systHigh['630GeV']['B'] = gTmp.GetErrorYhigh(
    0) / 3  # measurement integrated in |y|<0.6
extrapSystLow['630GeV']['B'] = 0.
extrapSystHigh['630GeV']['B'] = 0.
inFileUA1.Close()

inFileCDF = TFile.Open('ptintcross/HEPData_bbbar_CDF_ppbar1dot96TeV_PRD_71_032001_2005.root')
gTmp = inFileCDF.Get('Table 5/Graph1D_y1')
gTmp.GetPoint(0, E, dsdy)
sigma['1dot96TeV'], stat['1dot96TeV'], systLow['1dot96TeV'], systHigh['1dot96TeV'], \
    extrapSystLow['1dot96TeV'], extrapSystHigh['1dot96TeV'] = ({} for _ in range(6))
sigma['1dot96TeV']['Jpsi'] = dsdy.value / 1.2 # measurement integrated in |y|<0.6
stat['1dot96TeV']['Jpsi'] = 0. # info not available
systLow['1dot96TeV']['Jpsi'] = gTmp.GetErrorYlow(0) / 1.2
systHigh['1dot96TeV']['Jpsi'] = gTmp.GetErrorYhigh(0) / 1.2
extrapSystLow['1dot96TeV']['Jpsi'] = 0.
extrapSystHigh['1dot96TeV']['Jpsi'] = 0.
inFileCDF.Close()

inFileJpsi = TFile.Open('ptintcross/HEPData_bbbar_ALICEpp7TeV_npJpsi_JHEP_1211_065_2012.root')
gTmp = inFileJpsi.Get('Table 6/Graph1D_y1')
gTmp.GetPoint(0, E, dsdy)
sigma['7TeV']['Jpsi'] = dsdy.value  # measurement integrated in |y|<0.6
stat['7TeV']['Jpsi'] = 0.  # info not available
systLow['7TeV']['Jpsi'] = gTmp.GetErrorYlow(0)
systHigh['7TeV']['Jpsi'] = gTmp.GetErrorYhigh(0)
extrapSystLow['7TeV']['Jpsi'] = 0.
extrapSystHigh['7TeV']['Jpsi'] = 0.
inFileJpsi.Close()

inFileD5TeV = TFile.Open(
    'ptintcross/HEPData_bbbar_ALICEpp5TeV_npD_JHEP_05_2021_220.root')
hTmp = inFileD5TeV.Get('Table 16/Hist1D_y1')
hTmpStat = inFileD5TeV.Get('Table 16/Hist1D_y1_e1')
hTmpSys1 = inFileD5TeV.Get('Table 16/Hist1D_y1_e2')
hTmpSys2 = inFileD5TeV.Get('Table 16/Hist1D_y1_e3')
hTmpSys3 = inFileD5TeV.Get('Table 16/Hist1D_y1_e4')
hTmpSys4 = inFileD5TeV.Get('Table 16/Hist1D_y1_e6')
hTmpSysExtrLow = inFileD5TeV.Get('Table 16/Hist1D_y1_e5minus')
hTmpSysExtrhigh = inFileD5TeV.Get('Table 16/Hist1D_y1_e5plus')
sigma['5dot02TeV']['av'] = hTmp.GetBinContent(
    1)  # measurement integrated in |y|<0.6
stat['5dot02TeV']['av'] = hTmpStat.GetBinContent(1)  # info not available
systLow['5dot02TeV']['av'] = np.sqrt(hTmpSys1.GetBinContent(
    1)**2 + hTmpSys2.GetBinContent(1)**2 + hTmpSys3.GetBinContent(1)**2 + hTmpSys4.GetBinContent(1)**2)
systHigh['5dot02TeV']['av'] = systLow['5dot02TeV']['av']
extrapSystLow['5dot02TeV']['av'] = -hTmpSysExtrLow.GetBinContent(1)
extrapSystHigh['5dot02TeV']['av'] = hTmpSysExtrhigh.GetBinContent(1)
inFileD5TeV.Close()

# compilation plot
gStat, gSyst, gSystExtrap = ({} for _ in range(3))
for meas in sigma['13TeV']:
    gStat[meas], gSyst[meas], gSystExtrap[meas] = (
        TGraphAsymmErrors(1) for _ in range(3))
    gStat[meas].SetPoint(0, sigma['13TeV'][meas], order[meas])
    gSyst[meas].SetPoint(0, sigma['13TeV'][meas], order[meas])
    gSystExtrap[meas].SetPoint(0, sigma['13TeV'][meas], order[meas])
    gStat[meas].SetPointError(0, stat['13TeV'][meas],
                              stat['13TeV'][meas], 0., 0.)
    gSyst[meas].SetPointError(
        0, systLow['13TeV'][meas], systHigh['13TeV'][meas], 0.12, 0.12)
    if extrapSystLow['13TeV'][meas] > 0:
        gSystExtrap[meas].SetPointError(0, extrapSystLow['13TeV'][meas], extrapSystHigh['13TeV'][meas],
                                        0.1, 0.1)
    else:
        gSystExtrap[meas] = None

    SetObjectStyle(gStat[meas], color=kBlack, fillstyle=0)
    SetObjectStyle(gSyst[meas], color=kBlack, fillstyle=0)
    if gSystExtrap[meas]:
        SetObjectStyle(gSystExtrap[meas], color=kGray+1)

# vs sqrts
gStatVsSqrts, gSystVsSqrts, gTotUncWoExtrapVsSqrts, gSystExtrapVsSqrts = (
    {} for _ in range(4))
for energy in sigma:
    gStatVsSqrts[energy], gSystVsSqrts[energy], gTotUncWoExtrapVsSqrts[energy], \
        gSystExtrapVsSqrts[energy] = ({} for _ in range(4))
    for meas in sigma[energy]:
        gStatVsSqrts[energy][meas], gSystVsSqrts[energy][meas], gTotUncWoExtrapVsSqrts[energy][meas], \
            gSystExtrapVsSqrts[energy][meas] = (
                TGraphAsymmErrors(1) for _ in range(4))
        gStatVsSqrts[energy][meas].SetPoint(
            0, energies[energy]*shift[meas], sigma[energy][meas])
        gSystVsSqrts[energy][meas].SetPoint(
            0, energies[energy]*shift[meas], sigma[energy][meas])
        gTotUncWoExtrapVsSqrts[energy][meas].SetPoint(
            0, energies[energy]*shift[meas], sigma[energy][meas])
        gSystExtrapVsSqrts[energy][meas].SetPoint(
            0, energies[energy]*shift[meas], sigma[energy][meas])
        gStatVsSqrts[energy][meas].SetPointError(
            0, 0., 0., stat[energy][meas], stat[energy][meas])
        gSystVsSqrts[energy][meas].SetPointError(0, energies[energy]*0.05, energies[energy]*0.05,
                                                 systLow[energy][meas], systHigh[energy][meas])
        gTotUncWoExtrapVsSqrts[energy][meas].SetPointError(0, energies[energy]*0.05, energies[energy]*0.05,
                                                           np.sqrt(
                                                               stat[energy][meas]**2 + systLow[energy][meas]**2),
                                                           np.sqrt(stat[energy][meas]**2 + systHigh[energy][meas]**2))

        if extrapSystLow[energy][meas] > 0:
            gSystExtrapVsSqrts[energy][meas].SetPointError(0, energies[energy]*0.04, energies[energy]*0.04,
                                                           extrapSystLow[energy][meas], extrapSystHigh[energy][meas])
        else:
            gSystExtrapVsSqrts[energy][meas] = None

        if meas not in ['av', 'B0'] and not (meas == 'Dzero' and energy == '13TeV'):
            SetObjectStyle(gStatVsSqrts[energy][meas], color=kBlack,
                           fillstyle=0, markerstyle=markerstyles[energy][meas])
            SetObjectStyle(gSystVsSqrts[energy][meas],
                           color=kBlack, fillstyle=0)
            SetObjectStyle(gTotUncWoExtrapVsSqrts[energy][meas], color=kBlack, fillstyle=0,
                           markerstyle=markerstyles[energy][meas], markersize=1.4)
            if gSystExtrapVsSqrts[energy][meas]:
                SetObjectStyle(gSystExtrapVsSqrts[energy][meas], color=kGray+1)
        elif meas == 'B0':
            SetObjectStyle(gStatVsSqrts[energy][meas], color=kRed,
                           fillstyle=0, markerstyle=markerstyles[energy][meas])
            SetObjectStyle(gSystVsSqrts[energy][meas],
                           color=kRed, fillstyle=0)
            SetObjectStyle(gTotUncWoExtrapVsSqrts[energy][meas], color=kRed, fillstyle=0,
                           markerstyle=markerstyles[energy][meas], markersize=1.4)
            if gSystExtrapVsSqrts[energy][meas]:
                SetObjectStyle(gSystExtrapVsSqrts[energy][meas], color=kRedMy)
        else:
            SetObjectStyle(gStatVsSqrts[energy][meas], color=kBlack, fillstyle=0,
                           markerstyle=kFullDiamond, markersize=1.4)
            SetObjectStyle(gSystVsSqrts[energy][meas], color=kBlack, fillstyle=0,
                           markerstyle=kFullDiamond, markersize=1.2)
            SetObjectStyle(gTotUncWoExtrapVsSqrts[energy][meas], color=kBlack, fillstyle=0,
                           markerstyle=kFullDiamond, markersize=1.4)
            if gSystExtrapVsSqrts[energy][meas]:
                SetObjectStyle(gSystExtrapVsSqrts[energy][meas], color=kGray+1, fillalpha=0.5)


# D mesons
with open('config_filenames.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
fileNames = {'Dzero': inputCfg['ptintcross']['bbbar']['Dzero'],
             'Dplus': inputCfg['ptintcross']['bbbar']['Dplus'],
             'Ds': inputCfg['ptintcross']['bbbar']['Ds'],
             'Lc': inputCfg['ptintcross']['bbbar']['Lc'],
             'DzeroLc': inputCfg['ptintcross']['bbbar']['DzeroLc'],
             'av': inputCfg['ptintcross']['bbbar']['average']}

for hadron in fileNames:
    inFile = TFile.Open(fileNames[hadron])
    hStatTmp = inFile.Get('hbbbarCrossSecStat')
    gSystDataTmp = inFile.Get('gbbbarCrossSecDataSyst')
    gSystExtrapTmp = inFile.Get('gbbbarCrossSecExtrapSyst')
    gSystExtrapRapTmp = inFile.Get('gbbbarCrossSecRapExtrapSyst')

    gStat[hadron], gSyst[hadron], gSystExtrap[hadron] = (
        TGraphAsymmErrors(1) for _ in range(3))
    gStat[hadron].SetPoint(0, hStatTmp.GetBinContent(1), order[hadron]-0.1)
    gSyst[hadron].SetPoint(0, hStatTmp.GetBinContent(1), order[hadron]-0.1)
    gSystExtrap[hadron].SetPoint(
        0, hStatTmp.GetBinContent(1), order[hadron]-0.1)
    gStat[hadron].SetPointError(0, hStatTmp.GetBinError(
        1), hStatTmp.GetBinError(1), 0., 0.)
    gSyst[hadron].SetPointError(0, gSystDataTmp.GetErrorYlow(
        0), gSystDataTmp.GetErrorYhigh(0), 0.12, 0.12)
    gSystExtrap[hadron].SetPointError(0,
                                      np.sqrt(gSystExtrapTmp.GetErrorYlow(
                                          0)**2 + gSystExtrapRapTmp.GetErrorYlow(0)**2),
                                      np.sqrt(gSystExtrapTmp.GetErrorYhigh(
                                          0)**2 + gSystExtrapRapTmp.GetErrorYhigh(0)**2),
                                      0.1, 0.1)
    if hadron == 'av':
        SetObjectStyle(gStat[hadron], color=kBlack, fillstyle=0)
        SetObjectStyle(gSyst[hadron], color=kBlack, fillstyle=0)
        SetObjectStyle(gSystExtrap[hadron], color=kGray+1, fillalpha=0.5)
    elif hadron == 'DzeroLc':
        SetObjectStyle(gStat[hadron], color=kRed, fillstyle=0)
        SetObjectStyle(gSyst[hadron], color=kRed, fillstyle=0)
        SetObjectStyle(gSystExtrap[hadron], color=kRed)
    else:
        SetObjectStyle(gStat[hadron], color=kBlack, fillstyle=0)
        SetObjectStyle(gSyst[hadron], color=kBlack, fillstyle=0)
        SetObjectStyle(gSystExtrap[hadron], color=kGray+1, fillalpha=0.5)

    if hadron in ['av', 'Dzero', 'Lc', 'DzeroLc']:
        gStatVsSqrts['13TeV'][hadron], gSystVsSqrts['13TeV'][hadron], gTotUncWoExtrapVsSqrts['13TeV'][hadron], \
            gSystExtrapVsSqrts['13TeV'][hadron] = (
                TGraphAsymmErrors(1) for _ in range(4))
        gStatVsSqrts['13TeV'][hadron].SetPoint(
            0, energies['13TeV'], hStatTmp.GetBinContent(1))
        gSystVsSqrts['13TeV'][hadron].SetPoint(
            0, energies['13TeV'], hStatTmp.GetBinContent(1))
        gTotUncWoExtrapVsSqrts['13TeV'][hadron].SetPoint(
            0, energies['13TeV'], hStatTmp.GetBinContent(1))
        gSystExtrapVsSqrts['13TeV'][hadron].SetPoint(
            0, energies['13TeV'], hStatTmp.GetBinContent(1))
        gStatVsSqrts['13TeV'][hadron].SetPointError(
            0, 0., 0., hStatTmp.GetBinError(1), hStatTmp.GetBinError(1))
        gSystVsSqrts['13TeV'][hadron].SetPointError(0, energies['13TeV']*0.05, energies['13TeV']*0.05,
                                                    gSystDataTmp.GetErrorYlow(0), gSystDataTmp.GetErrorYhigh(0))
        gTotUncWoExtrapVsSqrts['13TeV'][hadron].SetPointError(0, energies['13TeV']*0.05,
                                                              energies['13TeV']*0.05,
                                                              np.sqrt(hStatTmp.GetBinError(
                                                                  1)**2 + gSystDataTmp.GetErrorYlow(0)**2),
                                                              np.sqrt(hStatTmp.GetBinError(1)**2 + gSystDataTmp.GetErrorYhigh(0)**2))
        gSystExtrapVsSqrts['13TeV'][hadron].SetPointError(0, energies['13TeV']*0.04, energies['13TeV']*0.04,
                                                          np.sqrt(gSystExtrapTmp.GetErrorYlow(
                                                              0)**2 + gSystExtrapRapTmp.GetErrorYlow(0)**2),
                                                          np.sqrt(gSystExtrapTmp.GetErrorYhigh(0)**2 + gSystExtrapRapTmp.GetErrorYhigh(0)**2))
        if hadron in ['av', 'Dzero']:
            SetObjectStyle(gStatVsSqrts['13TeV'][hadron], color=kBlack, fillstyle=0,
                           markerstyle=kFullDiamond, markersize=1.4)
            SetObjectStyle(gSystVsSqrts['13TeV'][hadron], color=kBlack, fillstyle=0,
                           markerstyle=kFullDiamond, markersize=1.2)
            SetObjectStyle(gTotUncWoExtrapVsSqrts['13TeV'][hadron], color=kBlack, fillstyle=0,
                           markerstyle=kFullDiamond, markersize=1.4)
            SetObjectStyle(gSystExtrapVsSqrts['13TeV'][hadron], color=kGray+1, fillalpha=0.5)
        elif hadron == 'Lc':
            SetObjectStyle(gStatVsSqrts['13TeV'][hadron], color=kOrange+7, fillstyle=0,
                           markerstyle=kFullCircle, markersize=1.4)
            SetObjectStyle(gSystVsSqrts['13TeV'][hadron], color=kOrange+7, fillstyle=0,
                           markerstyle=kFullCircle, markersize=1.2)
            SetObjectStyle(gTotUncWoExtrapVsSqrts['13TeV'][hadron], color=kOrange+7, fillstyle=0,
                           markerstyle=kFullCircle, markersize=1.4)
            SetObjectStyle(
                gSystExtrapVsSqrts['13TeV'][hadron], color=kOrange+6)
        elif hadron == 'DzeroLc':
            SetObjectStyle(gStatVsSqrts['13TeV'][hadron], color=kBlack, fillstyle=0,
                           markerstyle=kFullCross, markersize=1.4)
            SetObjectStyle(gSystVsSqrts['13TeV'][hadron], color=kBlack, fillstyle=0,
                           markerstyle=kFullCross, markersize=1.2)
            SetObjectStyle(gTotUncWoExtrapVsSqrts['13TeV'][hadron], color=kBlack, fillstyle=0,
                           markerstyle=kFullCross, markersize=1.4)
            SetObjectStyle(
                gSystExtrapVsSqrts['13TeV'][hadron], color=kGray+1, markerstyle=kFullCross)

    if hadron == 'B0':
        gStatVsSqrts['13dot6TeV'][hadron], gSystVsSqrts['13dot6TeV'][hadron], gTotUncWoExtrapVsSqrts['13dot6TeV'][hadron], \
            gSystExtrapVsSqrts['13dot6TeV'][hadron] = (
                TGraphAsymmErrors(1) for _ in range(4))
        gStatVsSqrts['13dot6TeV'][hadron].SetPoint(
            0, energies['13dot6TeV'], hStatTmp.GetBinContent(1))
        gSystVsSqrts['13dot6TeV'][hadron].SetPoint(
            0, energies['13dot6TeV'], hStatTmp.GetBinContent(1))
        gTotUncWoExtrapVsSqrts['13dot6TeV'][hadron].SetPoint(
            0, energies['13dot6TeV'], hStatTmp.GetBinContent(1))
        gSystExtrapVsSqrts['13dot6TeV'][hadron].SetPoint(
            0, energies['13dot6TeV'], hStatTmp.GetBinContent(1))
        gStatVsSqrts['13dot6TeV'][hadron].SetPointError(
            0, 0., 0., hStatTmp.GetBinError(1), hStatTmp.GetBinError(1))
        gSystVsSqrts['13dot6TeV'][hadron].SetPointError(0, energies['13dot6TeV']*0.05, energies['13dot6TeV']*0.05,
                                                    gSystDataTmp.GetErrorYlow(0), gSystDataTmp.GetErrorYhigh(0))
        gTotUncWoExtrapVsSqrts['13dot6TeV'][hadron].SetPointError(0, energies['13dot6TeV']*0.05,
                                                              energies['13dot6TeV']*0.05,
                                                              np.sqrt(hStatTmp.GetBinError(
                                                                  1)**2 + gSystDataTmp.GetErrorYlow(0)**2),
                                                              np.sqrt(hStatTmp.GetBinError(1)**2 + gSystDataTmp.GetErrorYhigh(0)**2))
        gSystExtrapVsSqrts['13dot6TeV'][hadron].SetPointError(0, energies['13dot6TeV']*0.04, energies['13dot6TeV']*0.04,
                                                          np.sqrt(gSystExtrapTmp.GetErrorYlow(
                                                              0)**2 + gSystExtrapRapTmp.GetErrorYlow(0)**2),
                                                          np.sqrt(gSystExtrapTmp.GetErrorYhigh(0)**2 + gSystExtrapRapTmp.GetErrorYhigh(0)**2))
        
        SetObjectStyle(gStatVsSqrts['13dot6TeV'][hadron], color=kRed, fillstyle=0,
                        markerstyle=kFullCross, markersize=1.4)
        SetObjectStyle(gSystVsSqrts['13dot6TeV'][hadron], color=kRed, fillstyle=0,
                        markerstyle=kFullCross, markersize=1.2)
        SetObjectStyle(gTotUncWoExtrapVsSqrts['13dot6TeV'][hadron], color=kRed, fillstyle=0,
                        markerstyle=kFullCross, markersize=1.4)
        SetObjectStyle(
            gSystExtrapVsSqrts['13dot6TeV'][hadron], color=kRedMy, markerstyle=kFullCross)


# FONLL
fonll = {}
gFONLL, gFONLLvsSqrts, gFONLLvsSqrtsCent = (
    TGraphAsymmErrors(0) for _ in range(3))

iE = 0
for energy in energies:
    if energy == '68dot6GeV' or energy == '200GeV' or energy == '630GeV' or energy == '900GeV' or energy == '1dot96TeV' or energy == '2dot76TeV':
        continue
    _, fonll[energy] = ReadFONLL(inputCfg['models']['FONLL']['bbbar'][energy])

    gFONLLvsSqrts.SetPoint(iE, energies[energy], fonll[energy]['central']/1.e6)
    gFONLLvsSqrtsCent.SetPoint(
        iE, energies[energy], fonll[energy]['central']/1.e6)
    gFONLLvsSqrts.SetPointError(iE, 0., 0., (fonll[energy]['central']-fonll[energy]['min'])/1.e6,
                                (fonll[energy]['max']-fonll[energy]['central'])/1.e6)

    if energy == '13TeV':
        gFONLL.SetPoint(0, fonll[energy]['central'] /
                        1.e6, len(measToPlot13TeV)/2+0.5)
        gFONLL.SetPointError(0, (fonll[energy]['central']-fonll[energy]['min'])/1.e6,
                             (fonll[energy]['max'] -
                              fonll[energy]['central'])/1.e6,
                             len(measToPlot13TeV)/2, len(measToPlot13TeV)/2)
    iE += 1

SetObjectStyle(gFONLL, linecolor=colors[0], fillcolor=colors[0], fillstyle=1001, alpha=0.5)
SetObjectStyle(gFONLLvsSqrts, color=colors[0], alpha=0.5, fillstyle=1001)
SetObjectStyle(gFONLLvsSqrtsCent, color=colors[0], alpha=0.5)
lFONLL = TLine(fonll['13TeV']['central']/1.e6, 0.5,
               fonll['13TeV']['central']/1.e6, len(measToPlot13TeV)+0.5)
lFONLL.SetLineColorAlpha(colors[0], 0.5)
lFONLL.SetLineWidth(2)

# NNLO
# pp 200 GeV: 0.750(12)  -20%  +21%
# pp 500 GeV: 3.226(29)  -20%  +19%
# pp 2.75 TeV: 20.52(17) -22%   +24%
# pp 5.02 TeV: 34.5(3)  -23% +28%
# pp 7 TeV: 44.69(24) -23% + 30%
# pp 13 TeV: 72.2(6)   -27% +31%
# ppbar 900 GeV:  6.41(17) -19% +12%
# ppbar 1.96 TeV: 15.2(3) -22% + 21%

nnlo = {}
gNNLO, gNNLOvsSqrts, gNNLOvsSqrtsCent = (
    TGraphAsymmErrors(0) for _ in range(3))
for energy in energies:
    nnlo[energy] = {}

nnlo['200GeV']['central'] = 0.75012
nnlo['200GeV']['min'] = 0.75012*(1-0.20)
nnlo['200GeV']['max'] = 0.75012*(1+0.21)

nnlo['500GeV']['central'] = 3.22629
nnlo['500GeV']['min'] = 3.22629*(1-0.20)
nnlo['500GeV']['max'] = 3.22629*(1+0.19)

nnlo['900GeV']['central'] = 6.4117
nnlo['900GeV']['min'] = 6.4117*(1-0.19)
nnlo['900GeV']['max'] = 6.4117*(1+0.12)

nnlo['1dot96TeV']['central'] = 15.23
nnlo['1dot96TeV']['min'] = 15.23*(1-0.22)
nnlo['1dot96TeV']['max'] = 15.23*(1+0.21)

nnlo['2dot76TeV']['central'] = 20.5217
nnlo['2dot76TeV']['min'] = 20.5217*(1-0.22)
nnlo['2dot76TeV']['max'] = 20.5217*(1+0.24)

nnlo['5dot02TeV']['central'] = 34.53
nnlo['5dot02TeV']['min'] = 34.53*(1-0.22)
nnlo['5dot02TeV']['max'] = 34.53*(1+0.21)

nnlo['7TeV']['central'] = 44.6924
nnlo['7TeV']['min'] = 44.6924*(1-0.23)
nnlo['7TeV']['max'] = 44.6924*(1+0.30)

nnlo['13TeV']['central'] = 72.26
nnlo['13TeV']['min'] = 72.26*(1-0.27)
nnlo['13TeV']['max'] = 72.26*(1+0.31)

nnlo['13dot6TeV']['central'] = 72.26
nnlo['13dot6TeV']['min'] = 72.26*(1-0.27)
nnlo['13dot6TeV']['max'] = 72.26*(1+0.31)

iE = 0
for energy in energies:
    if energy == '200GeV' or energy == '630GeV' or energy == '68dot6GeV' or energy == '5dot36TeV' or energy == '5dot5TeV' or energy == '8TeV' or energy == '8dot16TeV':
        continue

    gNNLOvsSqrts.SetPoint(iE, energies[energy], nnlo[energy]['central'])
    gNNLOvsSqrtsCent.SetPoint(iE, energies[energy], nnlo[energy]['central'])
    gNNLOvsSqrts.SetPointError(iE, 0., 0., (nnlo[energy]['central']-nnlo[energy]['min']),
                               (nnlo[energy]['max']-nnlo[energy]['central']))

    if energy == '13TeV':
        gNNLO.SetPoint(0, nnlo[energy]['central'], len(measToPlot13TeV)/2+0.5)
        gNNLO.SetPointError(0, (nnlo[energy]['central']-nnlo[energy]['min']),
                            (nnlo[energy]['max']-nnlo[energy]['central']),
                            len(measToPlot13TeV)/2, len(measToPlot13TeV)/2)
    iE += 1

SetObjectStyle(gNNLO, linecolor=colors[1], fillcolor=colors[1], alpha=0.5,
               fillstyle=1001, linestyle=9)
SetObjectStyle(gNNLOvsSqrts, fillcolor=colors[1], fillstyle=1001, fillalpha=0.5)
SetObjectStyle(gNNLOvsSqrtsCent, linecolor=colors[1], linestyle=9, alpha=0.5)

lNNLO = TLine(nnlo['13TeV']['central'], 0.5,
              nnlo['13TeV']['central'], len(measToPlot13TeV)+0.5)
lNNLO.SetLineColor(colors[1])
lNNLO.SetLineWidth(2)
lNNLO.SetLineStyle(9)

# draw results
SetGlobalStyle(padleftmargin=0.15, padbottommargin=0.14,
               labelsize=0.045, titlesizex=0.05)
cDiffMeas = TCanvas('cDiffMeas', '', 800, 600)
cSqrts = TCanvas('cSqrts', '', 800, 600)

latRef = TLatex()
latRef.SetTextFont(42)
latRef.SetTextColor(kBlack)
latRef.SetTextSize(0.035)

lat = TLatex()
lat.SetNDC()
lat.SetTextFont(42)
lat.SetTextColor(kBlack)
lat.SetTextSize(0.04)

latALICE = TLatex()
latALICE.SetNDC()
latALICE.SetTextFont(42)
latALICE.SetTextColor(kBlack)
latALICE.SetTextSize(0.05)

leg = TLegend(0.7, 0.6, 0.93, 0.8)
leg.SetTextSize(0.04)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetMargin(0.5)
leg.AddEntry(gFONLL, 'FONLL', 'lf')
leg.AddEntry(gNNLO, 'NNLO', 'lf')
leg.AddEntry(gSyst['Dzero'], 'data sys', 'f')
leg.AddEntry(gSystExtrap['Dzero'], 'extrap sys', 'f')

latALICEPreliminary = TLatex()
latALICEPreliminary.SetNDC()
latALICEPreliminary.SetTextFont(42)
latALICEPreliminary.SetTextColor(kBlack)
latALICEPreliminary.SetTextSize(0.05)


legVsSqrts = TLegend(0.15, 0.73, 0.45, 0.88)
legVsSqrts.SetTextSize(0.05)
legVsSqrts.SetFillStyle(0)
legVsSqrts.SetBorderSize(0)
legVsSqrts.SetMargin(0.1)
#legVsSqrts.AddEntry(gTotUncWoExtrapVsSqrts['200GeV']['HFe'],
#                    'PHENIX pp, |#it{y}|<1.5', 'p')  # lower[-0.06]{#scale[0.72]{PRL 103 (2009) 082002}}
legVsSqrts.AddEntry(gTotUncWoExtrapVsSqrts['630GeV']['B'],
                    'UA1 p#bar{p}, |#it{y}|<1.5', 'p')  # lower[-0.06]{#scale[0.72]{PLB 256 (1991) 121-128}}
legVsSqrts.AddEntry(gTotUncWoExtrapVsSqrts['1dot96TeV']['Jpsi'],
                    'CDF p#bar{p}, |#it{y}|<0.6', 'p')  # lower[-0.06]{#scale[0.72]{PRD 71 (2005) 032001}}

legVsSqrtsALICE = TLegend(0.55, 0.17, 0.93, 0.47)
legVsSqrtsALICE.SetTextSize(0.05)
legVsSqrtsALICE.SetFillStyle(0)
legVsSqrtsALICE.SetBorderSize(0)
legVsSqrtsALICE.SetMargin(0.1)
legVsSqrtsALICE.SetNColumns(1)
legVsSqrtsALICE.SetHeader('ALICE pp', 'R')
legVsSqrtsALICE.AddEntry(
    gTotUncWoExtrapVsSqrts['13dot6TeV']['B0'], 'b#rightarrow B^{0} |#it{y}|<0.5', 'p')
legVsSqrtsALICE.AddEntry(
    gTotUncWoExtrapVsSqrts['13TeV']['av'], 'b#rightarrow D |#it{y}|<0.5', 'p')
#legVsSqrtsALICE.AddEntry('', 'dielectron |#it{#eta}_{e}|<0.8', '')
legVsSqrtsALICE.AddEntry(
    gTotUncWoExtrapVsSqrts['13TeV']['DzeroLc'], 'b#rightarrow D + b#rightarrow#Lambda_{c}^{+} |#it{y}|<0.5', 'p')
#legVsSqrtsALICE.AddEntry('', '', '')
# legVsSqrtsALICE.AddEntry('', '#lower[-0.06]{#scale[0.72]{#sqrt{s} = 5.02, 13 TeV Preliminary}}', '')
legVsSqrtsALICE.AddEntry(
    gTotUncWoExtrapVsSqrts['13TeV']['Jpsi'], 'b#rightarrow J/#psi |#it{y}|<0.9', 'p')
#legVsSqrtsALICE.AddEntry('', '', '')
#legVsSqrtsALICE.AddEntry(gTotUncWoExtrapVsSqrts['7TeV']['HFe'], 'b#rightarrow e |#it{y}|<0.8', 'p')

# legVsSqrtsALICE.AddEntry('', '#lower[-0.06]{#scale[0.72]{#sqrt{s} = 7 TeV JHEP 1211 (2012) 065}}', '')
# legVsSqrtsALICE.AddEntry('', '  #lower[-0.06]{#scale[0.72]{#sqrt{s} = 5.02 TeV arXiv:2005.11995}}', '')
# legVsSqrtsALICE.AddEntry('', '  #lower[-0.06]{#scale[0.72]{#sqrt{s} = 7 TeV JHEP 1809 (2018) 064}}', '')
# legVsSqrtsALICE.AddEntry('', '#lower[-0.06]{#scale[0.72]{PLB 721 (2013) 13-23}}', '')
# legVsSqrtsALICE.AddEntry('', '  #lower[-0.06]{#scale[0.72]{#sqrt{s} = 13 TeV PLB 788 (2019) 505}}', '')

legDiel = TLegend(0.55, 0.15, 0.71, 0.4)
legDiel.SetTextSize(0.05)
legDiel.SetFillStyle(0)
legDiel.SetBorderSize(0)
legDiel.SetHeader('dielectron |#it{#eta}_{e}|<0.8')
legDiel.AddEntry(gTotUncWoExtrapVsSqrts['7TeV']['HFe'], 'b#rightarrow e |#it{y}|<0.8', 'p')
legDiel.AddEntry(
    gTotUncWoExtrapVsSqrts['13TeV']['eePy'], 'fit with PYTHIA6', 'p')
legDiel.AddEntry(
    gTotUncWoExtrapVsSqrts['13TeV']['eePO'], 'fit with POWHEG', 'p')

legFONLL = TLegend(0.15, 0.2, 0.61, 0.3)
legFONLL.SetTextSize(0.04)
legFONLL.SetFillStyle(0)
legFONLL.SetBorderSize(0)
legFONLL.AddEntry(gFONLL, 'FONLL', 'fl')
legFONLL.AddEntry(gNNLO, 'NNLO', 'fl')

cDiffMeas.cd()
hFrame = TH2F('hFrame', ';d#sigma_{b#bar{b}}/d#it{y}|_{#it{y}=0} (#mub);',
              100, 0., 200., len(measToPlot13TeV), 0.5, len(measToPlot13TeV)+0.5)
hFrame.GetYaxis().SetLabelSize(0.055)
for meas in measToPlot13TeV:
    hFrame.GetYaxis().SetBinLabel(order[meas], yTitles[meas])
cDiffMeas.SetLeftMargin(0.35)
hFrame.Draw()
gFONLL.Draw('2')
gNNLO.Draw('2')
lFONLL.Draw('same')
lNNLO.Draw('same')
for meas in measToPlot13TeV:
    if gSystExtrap[meas]:
        gSystExtrap[meas].Draw('2')
    gSyst[meas].Draw('2')
    gStat[meas].Draw('PZ')
    latRef.DrawLatex(refPosX[meas], order[meas]-0.5, references[meas])
hFrame.Draw('axissame')
latALICE.DrawLatex(0.76, 0.9, 'ALICE')
lat.DrawLatex(0.71, 0.85, 'pp, #sqrt{s} = 13 TeV')
leg.Draw()
cDiffMeas.Modified()
cDiffMeas.Update()

cDiffMeas.SaveAs('CrossSec_bbbar_compilation_pp13TeV.eps')
cDiffMeas.SaveAs('CrossSec_bbbar_compilation_pp13TeV.eps')
cDiffMeas.SaveAs('CrossSec_bbbar_compilation_pp13TeV.pdf')
cDiffMeas.SaveAs('CrossSec_bbbar_compilation_pp13TeV.pdf')

cSqrts.cd()
hFrameSqrts = cSqrts.DrawFrame(0.4, 0.5, 20., 1.5e2,
                               ';#sqrt{#it{s}} (TeV); d#sigma_{b#bar{b}}/d#it{y}|_{#it{y}=0} (#mub)')
hFrameSqrts.GetXaxis().SetMoreLogLabels()
hFrameSqrts.GetYaxis().SetTitleOffset(1)
# hFrameSqrts.GetXaxis().SetLabelSize(0.035)
cSqrts.SetLeftMargin(0.12)
cSqrts.SetLogy()
cSqrts.SetLogx()
gFONLLvsSqrts.Draw('3')
gNNLOvsSqrts.Draw('3')
gFONLLvsSqrtsCent.Draw('L')
gNNLOvsSqrtsCent.Draw('L')
for energy in measToPlotVsSqrt:
    for meas in measToPlotVsSqrt[energy]:
        if "ee" in meas or ("HFe" in meas and energy =="7TeV"):
            continue
        if gSystExtrapVsSqrts[energy][meas]:
            gSystExtrapVsSqrts[energy][meas].Draw('2')
        gTotUncWoExtrapVsSqrts[energy][meas].Draw('PZ')

legVsSqrts.Draw()
legVsSqrtsALICE.Draw()
#legDiel.Draw()
legFONLL.Draw()

latALICEPreliminary.DrawLatex(0.16, 0.89, 'ALICE Preliminary')
cSqrts.Modified()
cSqrts.Update()

cSqrts.SaveAs('CrossSec_bbbar_vs_sqrts.eps')
cSqrts.SaveAs('CrossSec_bbbar_vs_sqrts.pdf')

input('Press enter to exit')

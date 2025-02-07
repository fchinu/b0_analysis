import pandas as pd 
import ROOT

SIZE = 100000

file = ROOT.TFile.Open('simulations/b_to_dka/MDPi_511_from_DKa_decay_with_pt_cuts.root')
tree = file.Get('treeBtoDKa')

pt_mins = [1, 2, 4, 6, 8, 10, 14]
pt_maxs = [2, 4, 6, 8, 10, 14, 23.5]

df = pd.DataFrame(columns=['fM', 'fPt', 'ML_output'])
for pt_min, pt_max in zip(pt_mins, pt_maxs):
    tree.Draw('mBSmeared >> h_mb_smeared(760,4.9,5.66)', f'pTB>{pt_min} && pTB<{pt_max} && selected==1', 'goff')
    h_mb_smeared = ROOT.gDirectory.Get('h_mb_smeared')

    mass = []
    for i in range(SIZE):
        mass.append(h_mb_smeared.GetRandom())
    pt = [0.5 * (pt_min + pt_max)] * SIZE
    ml_score = [1] * SIZE
    b_mother = [511] * SIZE
    c_mother = [411+321] * SIZE
    df = pd.concat([df, pd.DataFrame({'fM': mass, 'fPt': pt,
                                    'ML_output': ml_score, 'fPdgCodeBeautyMother': b_mother, 'fPdgCodeCharmMother': c_mother})])

df.to_parquet("simulations/b_to_dka/dummy_for_template_with_pt_cuts.parquet")
    

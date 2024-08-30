def EvaluateEfficiencyFromHistos(hGen, hReco):
    """
    Evaluate the efficiency of a B0 particle.

    Args:
        hGen (TH1F): The histogram of the generated B0 particles.
        hReco (TH1F): The histogram of the reconstructed B0 particles.

    Returns:
        TH1F: The efficiency of the B0 particle.
    """
    
    hEff = hReco.Clone('hEff')
    hEff.Divide(hReco, hGen, 1., 1., 'B')
    return hEff

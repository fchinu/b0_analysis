import uproot
import ROOT
import numpy as np

def evaluate_efficiency_from_histos(hGen, hReco):
    """
    Evaluate the efficiency of a B0 particle.

    Args:
        hGen (TH1F): The histogram of the generated B0 particles.
        hReco (TH1F): The histogram of the reconstructed B0 particles.

    Returns:
        TH1F: The efficiency of the B0 particle.
    """
    
    hEff = hReco.Clone('h_efficiency')
    hEff.Divide(hReco, hGen, 1., 1., 'B')
    return hEff

def get_n_events_from_zorro(infile_names, zorro_folder, triggers_of_interest_names, h_collisions_path=None):
    """
    Get the total number of events from analysis output.

    Args:
        infile_names (str or list): The name(s) of the .root input file(s).
        zorro_folder (str): The name of the Zorro folder.
        triggers_of_interest_names (str or list): The name(s) of the triggers of interest.
        h_collisions_path (str): The path to the hCollisions histogram.

    Returns:
        float: The total number of events.
    """

    n_events = 0
    if not isinstance(infile_names, list):
        infile_names = [infile_names]

    if not isinstance(triggers_of_interest_names, list):
        triggers_of_interest_names = [triggers_of_interest_names]
        
    for infile_name in infile_names:
        f = uproot.open(infile_name)
        zorro = f[zorro_folder]

        for _, run_folder in zorro.items(recursive=False): # Loop over the run folders
            inspected_tvx = 0
            triggers_of_interest_skimming = 0
            analysed_triggers_of_interest = 0
            inspected_tvx += sum(run_folder['InspectedTVX'].values())
            triggers_of_interest_skimming += sum([run_folder['Selections'].values()[run_folder['Selections'].axis().labels().index(trigger_name)] for trigger_name in triggers_of_interest_names])
            analysed_triggers_of_interest += sum(run_folder['AnalysedTriggersOfInterest'].values())
            n_events += inspected_tvx*analysed_triggers_of_interest/triggers_of_interest_skimming
        
        if h_collisions_path is not None:
            h_collisions = f[h_collisions_path]
            z_vtx_eff = h_collisions.values()[h_collisions.axis().labels().index('PV #it{z}')]/h_collisions.values()[h_collisions.axis().labels().index('PV #it{z}') - 1]
            n_events = n_events*z_vtx_eff

    return n_events

def rebin_tgraph_asymm_errors(graph, new_bins):
    """
    Rebin a TGraphAsymmErrors object using a new binning scheme.
    The function checks if the new bins are compatible with the original bin edges.
    
    Parameters:
    - graph: TGraphAsymmErrors object to be rebinned.
    - new_bins: List or array of new bin edges.

    Returns:
    - new_graph: A new TGraphAsymmErrors object with the rebinned data.

    Raises:
    - ValueError: If the new bin edges are not compatible with the original bin edges.
    """
    # Convert new_bins to a numpy array for easier manipulation
    new_bins = np.array(new_bins)
    num_new_bins = len(new_bins) - 1

    # Extract the original bin edges from the graph
    original_bins = np.array([graph.GetX()[i] - graph.GetErrorXlow(i) for i in range(graph.GetN())] +
                             [graph.GetX()[graph.GetN() - 1] + graph.GetErrorXhigh(graph.GetN() - 1)])
    original_bins = np.unique(np.round(original_bins, decimals=8))  # Ensure no rounding issues

    # Check if the new bins are a subset of the original bins
    if not np.isin(new_bins, original_bins).all():
        raise ValueError("The new bins are not compatible with the original bin edges of the TGraphAsymmErrors.")

    # Initialize arrays to hold the new points and uncertainties
    new_x = np.zeros(num_new_bins)
    new_y = np.zeros(num_new_bins)
    new_x_err_low = np.zeros(num_new_bins)
    new_x_err_up = np.zeros(num_new_bins)
    new_y_err_low = np.zeros(num_new_bins)
    new_y_err_up = np.zeros(num_new_bins)

    # Loop over the new bins and calculate the rebinned values
    for i in range(num_new_bins):
        bin_mask = (graph.GetX() >= new_bins[i]) & (graph.GetX() < new_bins[i + 1])

        # Extract the relevant points
        x_vals = np.array([graph.GetX()[j] for j in range(graph.GetN()) if bin_mask[j]])
        y_vals = np.array([graph.GetY()[j] for j in range(graph.GetN()) if bin_mask[j]])
        x_err_low_vals = np.array([graph.GetErrorXlow(j) for j in range(graph.GetN()) if bin_mask[j]])
        x_err_up_vals = np.array([graph.GetErrorXhigh(j) for j in range(graph.GetN()) if bin_mask[j]])
        y_err_low_vals = np.array([graph.GetErrorYlow(j) for j in range(graph.GetN()) if bin_mask[j]])
        y_err_up_vals = np.array([graph.GetErrorYhigh(j) for j in range(graph.GetN()) if bin_mask[j]])

        # Compute the new x value as the weighted mean
        new_x[i] = np.average(x_vals)
        new_y[i] = np.sum(y_vals)/(new_bins[i + 1] - new_bins[i])

        # Propagate uncertainties
        new_y_err_low[i] = np.sqrt(np.sum(y_err_low_vals**2))/(new_bins[i + 1] - new_bins[i])
        new_y_err_up[i] = np.sqrt(np.sum(y_err_up_vals**2))/(new_bins[i + 1] - new_bins[i])

        # Assign the bin center and half-width as the new x error
        new_x_err_low[i] = new_x[i] - new_bins[i]
        new_x_err_up[i] = new_bins[i + 1] - new_x[i]

    # Create a new TGraphAsymmErrors with the rebinned data
    new_graph = ROOT.TGraphAsymmErrors(num_new_bins)
    for i in range(num_new_bins):
        new_graph.SetPoint(i, new_x[i], new_y[i])
        new_graph.SetPointError(i, new_x_err_low[i], new_x_err_up[i], new_y_err_low[i], new_y_err_up[i])

    return new_graph

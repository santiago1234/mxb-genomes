"""
Module with usefull functions
"""
import pandas as pd


def load_1tgp_metada(path_to_panel, filter_subpops=True):
    """
    Loads the 1TGP genomes sample information
    Args:
        path_to_panel: str, relative path to the file:
            resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel
        filter_subpops: bool, if yes subset the populations
            CHB|YRI|IBS|GBR|MXL|PEL|CLM|PUR
    Returns:
        pd.DataFrame, metadata
    """
    panel = pd.read_table(path_to_panel)
    panel = panel.iloc[:, :4]
    subpops = ['CHB', 'YRI', 'IBS', 'GBR', 'MXL', 'PEL', 'CLM', 'PUR']
    names = {
        'pop': 'Subpopulation',
        'super_pop': 'Superpopulation',
        'sample': 'Samplename'
    }
    panel.rename(names, axis=1, inplace=True)
    if filter_subpops:
        return panel[panel.Subpopulation.isin(subpops)]
    else:
        return panel

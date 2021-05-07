"""
Module with usefull functions
"""
import pandas as pd
import os


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
        'sample': 'Samplename',
        'gender': 'Gender'
    }
    panel.rename(names, axis=1, inplace=True)
    if filter_subpops:
        return panel[panel.Subpopulation.isin(subpops)]
    else:
        return panel


def load_mxb(path_to_mxb, full_info=False):
    """
    loads the MXB samples
    path_to_mxb: str, path to the file resources/genomes-metadata/50Genomes_info.txt
    """
    panel = pd.read_table(path_to_mxb)
    # add _ so the Samplename matches the samples
    # names in the vcf files.
    panel['Samplename'] = panel.Fam_ID.str.replace('MXB', 'MXB_')
    panel['Subpopulation'] = 'MXB'
    panel['Superpopulation'] = 'MXB'
    if full_info:
        return panel
    else:
        return panel.loc[:, ['Samplename', 'Subpopulation', 'Superpopulation', 'Gender']]


def load_populations_info(rel_path_to_project_root):
    """
    This function loads the population information for the
    populations that we are working with in this proyect.
    Args:
        rel_path_to_project_root: str, Relative path to
         the proyect root directory
    """
    rootp_to_1tgp = "resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel"
    rootp_to_mxb = "resources/genomes-metadata/50Genomes_info.txt"

    onetgp = load_1tgp_metada(os.path.join(rel_path_to_project_root, rootp_to_1tgp))
    mxb = load_mxb(os.path.join(rel_path_to_project_root, rootp_to_mxb))
    # change how geneder is encoded, so it is consitent
    mxb['Gender'] = mxb['Gender'].map({'F': 'female', 'M': 'male'})
    return pd.concat([onetgp, mxb], ignore_index=True)

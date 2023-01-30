"""
Compute the expected SFS under
the demographic model and get
the observer for the given variant
categorie.
"""

import gzip
import pickle
import click

import numpy as np
import pandas as pd
import moments
import demes


data = '../../data/220113-ConstructBoostrapedDatasets/data/whole-genome/spectrum-cat_intronic.pkl.gz'
mdl = '../220124-InfereModels/results/best-guest-NAT-EXPANSION-intronic.yml'

def load_data_and_model(data, mdl):
    '''

    Loads data and model

    Args:
        data: (str) path to file with joint SFS
        mdl: (str) path to inferred demographic model
    '''

    with gzip.open(data, "rb") as f:
        data = pickle.load(f)

    mdl = demes.load(mdl)
    data = data.project([50]*4)

    return data, mdl


def spectrum_fold_to_array(sf):
    '''Folds the expectrum and get a numpy array'''
    sf_folded = sf.fold()
    sf_folded = sf_folded[~sf_folded.mask].data
    return sf_folded


def to_frame(sfs_folded, pop, SFSfrom):
    '''
    Put the folded SFS in a pandas frame, with metadata info
    Args:
        sfs_folded: Spectrum
        pop: pop id
        SFSfrom: from model (expected) or data?
    '''
    minor_alle_f = list(range(1, len(sfs_folded)+1))

    d = {
        'Population': pop,
        'SF_from': SFSfrom,
        'Frequency': sfs_folded,
        'Minor_allel_freq': minor_alle_f
    }
    return pd.DataFrame(d)


def get_observed_and_expected_sfs(data, pop):
    index = data.pop_ids.index(pop)
    indices = list(range(4))
    indices.pop(index)
    sf_data = data.marginalize(indices)
    
    ## expected SFS under the model
    sf_expected = moments.Spectrum.from_demes(
        mdl,
        sampled_demes=[pop],
        sample_sizes=[50]
    )

    sf_expected = moments.Inference.optimally_scaled_sfs(sf_expected, sf_data)
    
    ## fold and put results in a data frame
    s_data = spectrum_fold_to_array(sf_data)
    s_expected = spectrum_fold_to_array(sf_expected)

    s_data = to_frame(s_data, pop, 'Data')
    s_expected = to_frame(s_expected, pop, 'Expected')

    return pd.concat([s_data, s_expected])


@click.command()
@click.option('--data', required=True)
@click.option('--mdl', required=True)
@click.option('--output', required=True)
@click.option('--category', required=True)
def main(data: str, mdl: str, output: str, category: str):
    """
    Computes Observed and Expected SFS.
    The Expected SFS is according to the model.

    Args:
        data: Joint SFS for the whole genome.
        mdl: yaml infferred model
        output: output filename for table
        category: Which variant category?
    """
    click.echo('loading data and model ...')
    data, mdl = load_data_and_model(data, mdl)
    results = pd.concat([get_observed_and_expected_sfs(data, x) for x in data.pop_ids])
    results['category'] = category
    results.to_csv(output, index=False)


if __name__ == '__main__':
    main()

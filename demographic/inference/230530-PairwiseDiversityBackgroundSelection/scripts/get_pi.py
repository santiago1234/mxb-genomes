"""
From Aaron:
    - Can you compute pairwise diversity (pi) for the different background selection partitions for intergenic and intronic data?
    - This can be done in moments using fs.pi()
    - Then by dividing that output by the uL from each subset of data, we get the relative differences in pairwise diversity across datasets
    - I want to see if the relative differences closely match the relative differences of $N_e$ inferred in the different datasets
"""
import gzip
import pickle

import pandas as pd
import moments

# Constant file path and project size for all functions
BASE_PATH = '../../data/230426-AnnotateChunksWithBvalues/data/whole-genome-quartile/q{}'
PROJECT_TO_SIZE = [40, 40, 40, 40]
POPULATIONS = ['YRI', 'IBS', 'CHB', 'MXB']


def load_ml_file(quartile, mut_type):
    """
    Load and process mL file.
    
    Parameters:
    quartile: the quartile of the data.
    mut_type: the mutation type, either 'intronic' or 'intergenic'.
    
    Returns:
    mL: the processed content of the mL file as a float.
    """
    file_name = 'mL_introns.txt' if mut_type == 'intronic' else 'mL_intergenic.txt'
    file_path = f'{BASE_PATH.format(quartile)}/{file_name}'

    with open(file_path, 'r') as file:
        mL = float(file.readline().replace('mL:', '').strip())

    return mL


def load_spectrum(quartile, mut_type):
    """
    Load the spectrum project and fold.
    
    Parameters:
    quartile: the quartile of the data.
    mut_type: the mutation type, either 'intronic' or 'intergenic'.
    """
    file_path = f'{BASE_PATH.format(quartile)}/spectrum-cat_{mut_type}.pkl.gz'

    print('Loading data...')
    with gzip.open(file_path, "rb") as file:
        spectrum = pickle.load(file)

    spectrum = spectrum.project(PROJECT_TO_SIZE)
    spectrum = spectrum.fold()

    return spectrum


def marginalize(spectrum, pop):
    """
    Get the marginal spectrum for a given population.
    
    Parameters:
    spectrum: The spectrum to be marginalized.
    pop: The population to keep in the marginalized spectrum.
    
    Returns:
    sf_pop: the marginalized spectrum.
    """
    # Get the indices of all populations except 'pop'
    indexes = [i for i, apop in enumerate(spectrum.pop_ids) if apop != pop]

    # Marginalize the spectrum
    sf_pop = spectrum.marginalize(indexes)

    return sf_pop


def get_pi(quartile, mut_type):
    """
    Get pi values for each population and return the data as a pandas DataFrame.
    
    Parameters:
    quartile: the quartile of the data.
    mut_type: the mutation type, either 'intronic' or 'intergenic'.
    
    Returns:
    data: A pandas DataFrame containing the data.
    """
    ml = load_ml_file(quartile, mut_type)
    sfs = load_spectrum(quartile, mut_type)

    # Get pi for each population
    pis = [marginalize(sfs, apop).pi() for apop in POPULATIONS]

    # Put data in a pandas DataFrame
    data = {
        'pop': POPULATIONS,
        'quartile': quartile,
        'mut_type': mut_type,
        'pi': pis,
        'ml': ml,
    }
    return pd.DataFrame(data)


results = [get_pi(quartile, mut_type) for quartile in range(1, 5) for mut_type in ['intronic', 'intergenic']]
results = pd.concat(results)
results.to_csv('results/pi.csv', index=False)

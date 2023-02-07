"""
The goal of this script is to compute the proportion of EAS
ancestry from the LAI results.
"""

from pathlib import Path
import sys

import pandas as pd

# Load tracts and some helper code I wrote
sys.path.append('../src/tracts-python3/')
sys.path.append('../../demographic/inference/220215-Tracts-CI/')

from tractsmodels.utils import list_individuals_in_dir
import tracts


populations = ['CLM', 'PEL', 'MXL', 'PUR']
labels = ['EUR', 'NAT', 'AFR', 'EAS']
TRACTS_DIR = '../210514-ProcessXGMixOutForTracts/data/4-pops/tracts/'


def get_tracts_data(population):
    """Tracts data

    Get the tracts data for the given population.

    Args:
        population (str): Population name 
    """
    inter = "_anc"
    end = "_cM.bed"

    dir_to_tracts = Path('.') / TRACTS_DIR / population
    assert dir_to_tracts.exists()

    individuals = list_individuals_in_dir(dir_to_tracts)
    poptracts = tracts.population(names=individuals, fname=(
        str(dir_to_tracts) + '/', inter, end))

    return poptracts


def mean_ancestry_propotions(population):
    """

    Compute ancestries mean proportions

    Args:
        population (str): Population name 
    """
    poptracts = get_tracts_data(population)
    anc_p = poptracts.get_mean_ancestry_proportions(ancestries=labels)
    anc_p = list(anc_p)
    anc_p.append(population)
    return anc_p


ancprops = [mean_ancestry_propotions(x) for x in populations]
ancprops = pd.DataFrame(ancprops, columns=labels + ['cohort'])
ancprops.to_csv('results/lai-ancestry-props.csv')

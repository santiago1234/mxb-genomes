"""
Helper functions to run the simulations
"""

from collections import namedtuple
import pandas as pd
import numpy as np
import msprime


# test_data

_files = {
    'region': '../220506-SetSimulationDesign/test-data/region_region_23.bed',
    'coding_region': '../220506-SetSimulationDesign/test-data/region_exons_23.bed',
    'noncoding_region': '../220506-SetSimulationDesign/test-data/region_intronANDinterg_23.bed',
    'ml_coding': '../220506-SetSimulationDesign/test-data/region_mlcoding_23.csv',
    'ml_noncoding': '../220506-SetSimulationDesign/test-data/region_mlnoncoding_23.txt'
}


_fields = [
        'chromosome',
        'start',
        'end',
        'coding_intervals',
        'noncoding_intervals',
        'ml_noncoding',
        'ml_synonymous',
        'ml_missense',
        'ml_LOF'
    ]

SimData = namedtuple("SimData", _fields)

def load_region(region_files, path_to_genetic_maps):
    """
    region_files: dict mapping names to filenames. Should
    contain the following key-value-pairs
        'region': The sampled 1Mb genomic region.
        'coding_region': Coding intervals.
        'noncoding_region': Intronic and intergenic intervals
        'ml_coding': Scaled mutation rates for missense, synonymous, and loss of function variants.
        'ml_noncoding': Scaled mutation rates for non-coding variants (intronic and intergenic).

    Returns:
        SimData
    """

    # Load the region for the simulation
    chromosome, start, end = np.loadtxt(region_files['region'], dtype=np.int64)

    # Load intervals where we will put the mutations
    coding_intervals = pd.read_csv(region_files['coding_region'], sep='\t', names=[
                                   'chro', 'start', 'end'])
    noncoding_intervals = pd.read_csv(
        region_files['noncoding_region'], sep='\t', names=['chro', 'start', 'end'])

    # Mls
    ml_coding = pd.read_csv(region_files['ml_coding'])
    ml_non_coding = np.loadtxt(
        region_files['ml_noncoding'], dtype=np.object0)[1]
    ml_non_coding = float(ml_non_coding)

    # Create a dict mapping class name to value

    ml = {}
    for x, y in ml_coding.iterrows():
        ml[y.Q] = y.mL

    ml['noncoding'] = ml_non_coding

    sim_data = SimData(
        chromosome=chromosome,
        start=start,
        end=end,
        coding_intervals=coding_intervals,
        noncoding_intervals=noncoding_intervals,
        ml_noncoding=ml['noncoding'],
        ml_missense=ml['missense'],
        ml_synonymous=ml['synonymous'],
        ml_LOF=ml['LOF']
    )

    # TODO: I need to start the region at 0

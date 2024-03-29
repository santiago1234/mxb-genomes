"""
Helper functions to run the simulations.

Here I define a class SimData to keep
the data needed to run one simulation.
"""

from collections import namedtuple
import os

import pandas as pd
import numpy as np
import msprime

from .simulation import load_sim_as_ts

# DFE
# These DFEs were inferred by Aaron, in this paper:
# https://academic.oup.com/genetics/advance-article/doi/10.1093/genetics/iyac097/6613932?login=true

_DFE_base = namedtuple('DFE', 'varian_class, shape, scale, h, Ne')


class DFE(_DFE_base):
    """
    Distribution of Fitnest Effects
    """
    __slots__ = ()

    def mean(self):
        """
        Get the mean of the distribution.
        """
        return (self.scale * self.shape) / (2 * self.Ne)


DFE_missense = DFE(varian_class='missense',
                   shape=0.147,
                   scale=2117,
                   h=0.5,
                   Ne=12300)
DFE_lof = DFE(varian_class='LOF', shape=0.188, scale=121419, h=0.5, Ne=12300)

# test_data

_files = {
    'region':
    '../220506-SetSimulationDesign/test-data/region_region_23.bed',
    'coding_region':
    '../220506-SetSimulationDesign/test-data/region_exons_23.bed',
    'noncoding_region':
    '../220506-SetSimulationDesign/test-data/region_intronANDinterg_23.bed',
    'ml_coding':
    '../220506-SetSimulationDesign/test-data/region_mlcoding_23.csv',
    'ml_noncoding':
    '../220506-SetSimulationDesign/test-data/region_mlnoncoding_23.txt'
}

_fields = [
    'chromosome', 'start', 'end', 'coding_intervals', 'noncoding_intervals',
    'rmap', 'ml_noncoding', 'ml_synonymous', 'ml_missense', 'ml_LOF'
]

BaseSimData = namedtuple("SimData", _fields)


class SimData(BaseSimData):
    """
    Class to hold simulation data.
    This object contains the data for running the simulation:
        - chromosome: The region falls in a chromosome
        - start: Start genomic position
        - end: End genomic position
        - coding_intervals: pd.DataFrame with intervals, the position is relative to start
        - noncoding_intervals: pd.DataFrame with intervals, the position is relative to start
        - rmap: the recombination map in the region
        - ml_noncoding: Scaled mutation rate (intergenic and intronic intervals)
        - ml_synonymous: Scaled mutation rate
        - ml_missense: Scaled mutation rate
        - ml_LOF: Scaled mutation rate
    """
    __slots__ = ()

    def __repr__(self):
        return f"Region: {self.chromosome}, start: {self.start}, end: {self.end}"

    def L(self, coding=True):
        """compute the length of the intervals
        either coding = True or non coding (coding=False) 
        """
        if coding:
            return (self.coding_intervals.end -
                    self.coding_intervals.start).sum()
        else:
            return (self.noncoding_intervals.end -
                    self.noncoding_intervals.start).sum()

    @property
    def m_noncoding(self):
        return self.ml_noncoding / self.L(coding=False)

    @property
    def m_synonymous(self):
        return self.ml_synonymous / self.L(coding=True)

    @property
    def m_missense(self):
        return self.ml_missense / self.L(coding=True)

    @property
    def m_LOF(self):
        return self.ml_LOF / self.L(coding=True)


def recombination_map(path_to_genetic_maps, chromosome, start, end):
    """
    Args:
        path_to_genetic_maps: path to mxb-genomes/resources/genetic-maps/
    """
    map_file_name = f'chr{chromosome}.b38.gmap'
    map_file_name = os.path.join(path_to_genetic_maps, map_file_name)
    rmap = msprime.RateMap.read_hapmap(map_file_name,
                                       position_col=0,
                                       map_col=2)
    # we can take a slice from the map to get the coordinates in the sampled region
    # with set trim=True
    rmap = rmap.slice(left=start, right=end, trim=True)
    return rmap


def load_sim_files(sample_id, path_to_samples):
    """
    This function relies on the names that I have created.
    Each sample file has a unique sample id.
    Returns:

    contain the following key-value-pairs
        'region': The sampled 1Mb genomic region.
        'coding_region': Coding intervals.
        'noncoding_region': Intronic and intergenic intervals
        'ml_coding': Scaled mutation rates for missense, synonymous, and loss of function variants.
        'ml_noncoding': Scaled mutation rates for non-coding variants (intronic and intergenic).
    """
    region_files = {
        'region': f'region_region_{sample_id}.bed',
        'coding_region': f'region_exons_{sample_id}.bed',
        'noncoding_region': f'region_intronANDinterg_{sample_id}.bed',
        'ml_coding': f'region_mlcoding_{sample_id}.csv',
        'ml_noncoding': f'region_mlnoncoding_{sample_id}.txt'
    }

    return {
        x: os.path.join(path_to_samples, region_files[x])
        for x in region_files
    }


def simuldata(path_to_samples, sample_id, path_to_genetic_maps):
    """
    Args:
        path_to_samples: relative path to directory containing the sample regions.
            e.g. the path to  mxb-genomes/demographic/data/220404-SimulationData/data/samples
        sample_id: The sample id (a number)
        path_to_genetic_maps: the relative path to  mxb-genomes/resources/genetic-maps/
    Returns:
        SimData
    """

    # Load the region for the simulation

    region_files = load_sim_files(sample_id, path_to_samples)
    chromosome, start, end = np.loadtxt(region_files['region'], dtype=np.int64)

    # Load intervals where we will put the mutations
    coding_intervals = pd.read_csv(region_files['coding_region'],
                                   sep='\t',
                                   names=['chro', 'start', 'end'])
    noncoding_intervals = pd.read_csv(region_files['noncoding_region'],
                                      sep='\t',
                                      names=['chro', 'start', 'end'])

    # Position of the intervals is relative to the start genomic position
    coding_intervals['start'] -= start
    coding_intervals['end'] -= start
    noncoding_intervals['start'] -= start
    noncoding_intervals['end'] -= start

    # Mls
    ml_coding = pd.read_csv(region_files['ml_coding'])
    ml_non_coding = np.loadtxt(region_files['ml_noncoding'],
                               dtype=np.object0)[1]
    ml_non_coding = float(ml_non_coding)

    # Create a dict mapping class name to value

    ml = {}
    for x, y in ml_coding.iterrows():
        ml[y.Q] = y.mL

    ml['noncoding'] = ml_non_coding

    # recombination map
    rmap = recombination_map(path_to_genetic_maps, chromosome, start, end)

    sim_data = SimData(chromosome=chromosome,
                       start=start,
                       end=end,
                       coding_intervals=coding_intervals,
                       noncoding_intervals=noncoding_intervals,
                       rmap=rmap,
                       ml_noncoding=ml['noncoding'],
                       ml_missense=ml['missense'],
                       ml_synonymous=ml['synonymous'],
                       ml_LOF=ml['LOF'])

    return sim_data


def relpath_to_datafiles(path_to_root: str):
    """
    This is a helper function to load
    data files for simulation and simulation
    output analysis
    Args:
        path_to_root: relative path to project root dir /mxb-genomes
    """
    d_files = {
        'genetic_maps':
        'resources/genetic-maps/',
        'genome_masks':
        'resources/genome-masks/20160622.allChr.mask.bed',
        'sim_chunks_data':
        'demographic/data/220404-SimulationData/data/samples/',
        'sim_out':
        'demographic/simulations/220728-Simulation-DFE-Demography/results/simulations/',
        'demes_graph':
        'demographic/simulations/220422-fwdpy11-initial-test/ADMIXTURE-MXL.yml',
    }

    d_files = {
        fname: os.path.join(path_to_root, d_files[fname])
        for fname in d_files
    }

    return d_files


def load_simulation_output(path_to_root: str, sim_id: int):
    """
    Helper function to load simulation data,
        this function assumes that files for simulation
        input and outout have already been generated.
    Returns:
        tuple: (simdat, sim_ts, d_files)
            simdat: SimData object used to generate the current simulation
            sim_ts: Resulting tree sequence from simulation
            d_files: see relpath_to_datafiles function
    """

    d_files = relpath_to_datafiles(path_to_root)

    for f_name, f_path in d_files.items():
        assert os.path.exists(f_path), f'file/dir: {f_path} not found'

    sim_file = f"{d_files['sim_out']}sim-{sim_id}-pop.bin"
    assert os.path.exists(sim_file), f'file not found: {sim_file}'
    print(f'loading data for simulation: {sim_id} ...')

    sim_ts = load_sim_as_ts(sim_file, d_files['demes_graph'])

    simdat = simuldata(d_files['sim_chunks_data'], sim_id,
                       d_files['genetic_maps'])

    return simdat, sim_ts, d_files

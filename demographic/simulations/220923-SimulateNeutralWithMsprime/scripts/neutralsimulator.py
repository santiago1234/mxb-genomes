"""
Run a neutral simulation, on the simulation data, using msprime
usage:
    python scripts/neutralsimulator.py <SIMID> <DEMES_GRAPH> <OUTFOLDER>
"""
import os
import sys

import demes
import msprime

sys.path.append('../')

from simutils import utils, simulation

POPs = ['YRI', 'IBS', 'CHB', 'MXB', 'MXL']
N_IND = 50


def define_samples():
    """
    Define the samples to simulate.
    NOTE:
        In the tree sequence the first
        N_IND will correspondo to population POPs[0],
        the 2nd N_IND to population POPs[1], and so on
    """
    return [msprime.SampleSet(N_IND, population=pop) for pop in POPs]


def simulate_neutral(simdata, graph):
    """
    Runs a neutral simulation using msprime
    Args:
        simdata: simulation data object
        graph: a demes graph, the demography
    Returns:
        ts_noncoding, ts_synonymous
    """
    print('setting simulation')
    demography = msprime.Demography.from_demes(graph)
    evolution_model = [
        msprime.DiscreteTimeWrightFisher(duration=20),
        msprime.StandardCoalescent()
    ]
    samples = define_samples()

    print(
        f'runing simulation at: chr{simdata.chromosome} {simdata.start:,} - {simdata.end:,}'
    )
    ts_ancestry = msprime.sim_ancestry(samples=samples,
                                       demography=demography,
                                       recombination_rate=simdata.rmap,
                                       model=evolution_model,
                                       ploidy=2,
                                       random_seed=42)

    # get the mutation rate maps
    nr_maps = simulation.neutral_rate_maps(simdata)

    ts_noncoding = msprime.sim_mutations(ts_ancestry,
                                         rate=nr_maps['noncoding'],
                                         random_seed=42)

    ts_synonymous = msprime.sim_mutations(ts_ancestry,
                                          rate=nr_maps['synonymous'],
                                          random_seed=42)

    print(f'No of noncoding mutations: {ts_noncoding.num_mutations}')
    print(f'No of synonymous mutations: {ts_synonymous.num_mutations}')
    return ts_noncoding, ts_synonymous


def main(simid, graph, outfolder):
    """
    Args:
        simid: simulation region ID
        graph: path to demes graph *.yml
        outfolder: output folder
    """

    path_to_samples = '../../data/220404-SimulationData/data/samples/'
    path_to_genetic_maps = '../../../resources/genetic-maps/'
    graph = demes.load(graph)

    simdata = utils.simuldata(path_to_samples, simid, path_to_genetic_maps)

    outfiles = [
        f'msprime_sim_{simid}-cat_{x}.ts' for x in ('noncoding', 'synonymous')
    ]
    outfiles = [os.path.join(outfolder, x) for x in outfiles]

    ts_noncod, ts_syn = simulate_neutral(simdata, graph)

    ts_noncod.dump(outfiles[0])
    ts_syn.dump(outfiles[1])


if __name__ == "__main__":
    main(*sys.argv[1:])

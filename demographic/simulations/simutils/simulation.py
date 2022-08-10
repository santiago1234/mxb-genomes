"""
Helpfull code to set the simulation
"""
import pybedtools
import pandas as pd
import numpy as np
import msprime
import moments
import fwdpy11
import tskit
import demes
# *************** ========== ***************
# NEUTRAL SIMULATION: RATE MAPS
# *************** ========== ***************


def _make_rate_map_table(region, intervals_with_positve_rate, rate):
    """
    The inpute intervals will have the given rate,
    the complement intervals (to cover the whole region) will have 0 rate

    Args:
        regions: chrm,start, end: the genomic region bedtools
        intervals_with_positve_rate: bedtools intervals within the region
    Return:
        pd.DataFrame
    """

    # get the complement regions that will have zero rate
    intervals_zero_rate = region.subtract(intervals_with_positve_rate)

    # put in data frame
    intervals_with_positve_rate = intervals_with_positve_rate.to_dataframe()
    intervals_zero_rate = intervals_zero_rate.to_dataframe()

    # add rates accordingly
    intervals_with_positve_rate['rate'] = rate
    intervals_zero_rate['rate'] = 0
    intervals = pd.concat([intervals_with_positve_rate, intervals_zero_rate])

    return intervals.sort_values('start')


def _to_msprime_rate_map(rmt):
    """
    Args:
        rmt: rate map table
    Return:
        msprime.RateMap
    """
    positions = rmt.start.to_list()
    # append the last position
    positions.append(rmt.end.to_list()[-1])
    return msprime.RateMap(position=positions, rate=rmt.rate.to_list())


def rate_map(simdat, synonymous=True):
    """
    Make an msprime rate map
    """

    if synonymous:
        u = simdat.m_synonymous
        regions_with_rate = simdat.coding_intervals
    else:  # intronic and intergenic regions
        u = simdat.m_noncoding
        regions_with_rate = simdat.noncoding_intervals

    # make Bedtool

    start = simdat.start - simdat.start
    end = simdat.end - simdat.start
    region = pybedtools.BedTool(
        f'{simdat.chromosome} {start} {end}', from_string=True)

    regions_with_rate = pybedtools.BedTool.from_dataframe(regions_with_rate)

    rmap_table = _make_rate_map_table(region, regions_with_rate, u)
    return _to_msprime_rate_map(rmap_table)


def neutral_rate_maps(simdat):
    """
    Get the rate map for neutral regions:
    coding and no
    """
    neutral_maps = dict()
    neutral_maps['synonymous'] = rate_map(simdat, synonymous=True)
    neutral_maps['noncoding'] = rate_map(simdat, synonymous=False)

    return neutral_maps


# *************** ========== ***************
# NEUTRAL SIMULATION: SIM NEUTRAL MUTATIONS
# *************** ========== ***************

def simulate_neutral_variation(ts, simdat):
    """
    Simulate neutral genetic variation (noncoding and synonymous)
    Args:
        ts: tree sequence from forwar in time simulation
        simdat: the simulation data used to generate ts
    Returns:
        (ts_noncoding, ts_synonymous): the tree sequences with the 
            neutral mutations only
    """
    # get a dict with the rate maps for nuetral categories
    # synonymous and noncoding
    nr_maps = neutral_rate_maps(simdat)

    # remove existing mutations (selected) from (ts)
    # get the list of all ids
    all_ids = [v.site.id for v in ts.variants()]
    ts_clear = ts.delete_sites(all_ids)

    ts_nocd = msprime.sim_mutations(
        tree_sequence=ts_clear, rate=nr_maps['noncoding'])
    ts_syn = msprime.sim_mutations(
        tree_sequence=ts_clear, rate=nr_maps['synonymous'])

    return ts_nocd, ts_syn


# *************** ========== ***************
# SIMULATION OUTPUT PROCESSING
# *************** ========== ***************

def load_sim_as_ts(popbin, graph):
    """
    Loads the simulation output as a tree sequence
    Args:
        popbin: str, path to simulation output
        graph: str, path to demes model graph
    Returns:
        treesequence
    """
    pop = fwdpy11.DiploidPopulation.load_from_file(popbin)
    # we load the model to include this info in the tree sequence
    graph = demes.load(graph)
    pop_md = {}
    for i, deme in enumerate(graph.demes):
        pop_md[i] = {'name': deme.name, "description": deme.description}
    ts = pop.dump_tables_to_tskit(
        demes_graph=graph, population_metadata=pop_md)
    return ts


def subsample_individuals_pop(ts, N):
    """
    Take a random sample of N individuals from each deme in the tree
    """
    # group indivuals by population
    inds_by_pop = dict()

    for ind in ts.individuals():
        current_deme = ind.metadata['deme']
        if current_deme in inds_by_pop.keys():
            inds_by_pop[current_deme].append(ind)
        else:
            inds_by_pop[current_deme] = [ind]

    # take the random sample, no replacement
    sampled_inds = []

    for x in inds_by_pop.keys():
        sample = np.random.choice(inds_by_pop[x], size=N, replace=False)
        sampled_inds.extend(sample)

    nodes_to_keep = []
    for ind in sampled_inds:
        nodes_to_keep.extend(ind.nodes)

    return ts.simplify(nodes_to_keep)

# *************** ========== ***************
# Compute Selected SFS
# *************** ========== ***************


mut_labels = {
    'neutral': 0,
    'missense': 1,
    'synonymous': 2,
    'LOF': 3,
}


def keep_selected_sites(ts, missense=True):
    """
    Keep only mutations that are missense or LOF
    The output treen will only contain missense or lof.
    """
    mut_type = 'missense' if missense else 'LOF'
    label = mut_labels[mut_type]
    sites = []

    for x in ts.mutations():
        if x.metadata['label'] != label:
            sites.append(x.site)

    return ts.delete_sites(sites)


def get_popname_to_popid_mapping(ts):
    """
    Returns a dict mapping population names
    to populations id.
    Args:
        ts: tree sequence containing population metadata
    Returns:
        dict
    """
    return {x.metadata['name']: x.id for x in ts.populations()}


def get_individuals_from_pops(ts, poplist):
    """
    Returns the list of nodes that belong to populations in poplist
        ts: tree sequence containing population metadata
        poplist: list of population names
    """
    pop_ids = get_popname_to_popid_mapping(ts)
    # the deme is the id if the pops we want
    demes = [pop_ids[x] for x in poplist]
    nodes_from_pops = []

    for ind in ts.individuals():
        if ind.metadata['deme'] in demes:
            nodes_from_pops.extend(ind.nodes)

    return nodes_from_pops


def get_single_pop_sfs(ts, pop):
    """
    Get the site frequency spectrum for the given pop
    """
    # subset ts to samples in pop
    nodes_from_pop = get_individuals_from_pops(ts, [pop])
    ts_pop = ts.simplify(nodes_from_pop)
    afs = ts_pop.allele_frequency_spectrum(
        polarised=True, span_normalise=False)
    return moments.Spectrum(afs)

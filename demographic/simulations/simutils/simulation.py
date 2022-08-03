"""
Helpfull code to set the simulation
"""
import pybedtools
import pandas as pd
import msprime

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
    # get a dict with the rate maps for nuetral categories
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

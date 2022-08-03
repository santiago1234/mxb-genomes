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

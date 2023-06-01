"""
Get the intervals
"""
from glob import glob
from os import path
import re

import pybedtools

N_intervals = 100

intervals = glob('../../data/220113-ConstructBoostrapedDatasets/data/chunks/chunk_*.bed')

interval = '../../data/220113-ConstructBoostrapedDatasets/data/chunks/chunk_107.bed'

def lines_in_interval(interval):
    """
    How many line are in the interval?
    """
    with open(interval) as f:
        return len(f.readlines())


# get intervals with only one line
intervals = [interval for interval in intervals if lines_in_interval(interval) == 1]


def interval_id(interval):
    """
    Each interval has a unique id
    this is in the filename, the number before the .bed extension
    """
    basename = path.basename(interval)
    return int(re.search('chunk_(\d+).bed', basename).group(1))


def load_noncoding_intervals(interval):
    """
    Load the non-coding intervals,
    associated with the interval.
    Intronic and intergenic intervals are merged.
    """
    intervalid = interval_id(interval)
    intergenic = f'../../data/220113-ConstructBoostrapedDatasets/data/mL-noncoding/regions/intergenic_chunk_{intervalid}.bed'
    intronic = f'../../data/220113-ConstructBoostrapedDatasets/data/mL-noncoding/regions/introns_chunk_{intervalid}.bed'
    intergenic = pybedtools.BedTool(intergenic)
    
    return intergenic.cat(intronic,postmerge=True)


for i, interval in enumerate(intervals):

    if i > N_intervals:
        break

    out_file = f'data/intervals/noncoding_intervals_{i}.bed'
    noncoding_intervals = load_noncoding_intervals(interval)
    noncoding_intervals.saveas(out_file)

"""
Get map slice for a given region (interval)

usage: python get_recmap.py <region_id>

Args:
    region_id: ID of the region to get the map for (a value from 1 to 350)
"""
import sys
import pandas as pd


def linear_interpolation(x1, y1, x2, y2, x):
    """
    Performs linear interpolation given two points (x1, y1) and (x2, y2),
    and the desired x-value (x).

    Returns the interpolated y-value.
    """
    # Calculate the slope (m)
    m = (y2 - y1) / (x2 - x1)

    # Calculate the y-intercept (b)
    b = y1 - m * x1

    # Perform linear interpolation
    y = m * x + b

    return y


def value_interpolation(rmap, position):
    """
    Get the interpolated cM value for a given position
    """

    int_l = rmap[rmap['pos'] < position].tail(1)
    int_r = rmap[rmap['pos'] > position].head(1)

    # Interpolate the cM value
    cM = linear_interpolation(
        int_l['pos'].to_list()[0],
        int_l['cM'].to_list()[0],
        int_r['pos'].to_list()[0],
        int_r['cM'].to_list()[0],
        position)
    return cM


def get_slice(rmap, start, end):
    """
    Get the slice of recombination map between start and end
    """

    rm_slice = rmap[(rmap['pos'] > start) & (rmap['pos'] < end)].copy()

    # Interpolation for end points
    chr_value = rmap['chr'].iloc[0]

    cM_start = value_interpolation(rmap, start)
    cM_end = value_interpolation(rmap, end)

    # Add the start and end points
    cM_start_df = pd.DataFrame({'pos': [start], 'cM': [cM_start], 'chr': [chr_value]})
    cM_end_df = pd.DataFrame({'pos': [end], 'cM': [cM_end], 'chr': [chr_value]})

    rm_slice = pd.concat([cM_start_df, rm_slice, cM_end_df])

    # Sort the slice
    rm_slice = rm_slice.sort_values(by='pos')

    # Subtract so it starts at 0
    rm_slice['cM'] = rm_slice['cM'] - cM_start

    # substract start from all positions
    rm_slice['pos'] = rm_slice['pos'] - start

    # set the chromosome to 1
    rm_slice['chr'] = 1

    return rm_slice


if __name__ == '__main__':
    region_id = sys.argv[1]
    region_path = f'../../simulations/220728-Simulation-DFE-Demography/results/simulations/sim-{region_id}-pop.bin'
    with open(region_path) as f:
        region = f.readline()

    chrom, start, end = map(int, region.strip().split('\t'))

    rmap_path = f'../../../resources/genetic-maps/chr{chrom}.b38.gmap'
    rmap = pd.read_csv(rmap_path, sep='\t')

    rslice = get_slice(rmap, start, end)
    outfile = 'data/recomb_map/sim_{region_id}-rmap.tsv'
    rslice = rslice[['pos', 'cM', 'chr']
    rslice.to_csv(outfile, sep='\t', index=False)


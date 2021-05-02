"""
Compute trackt length distribution
"""

import numpy as np
import pandas as pd


def is_ancestry_continuous(r_1, r_2):
    """
    Params:
        r_1: pd.DataFrame, a single row pandas frame.
            this row contains the coulummn Ancestry
        r_2: same as r_1
    Returns:
        lgl: True if r_1 and r_2 represent the same ancestry
    Note:
        I assume that r_2 follows from r_1. I need
        to write conditions to make sure that r_1 and r_2
        if merged together represent continuous chromosome
        fragments.
        Also if there is a gap of more than 10000, the
        fragments are not considered to be continous
    """
    max_gap = 10000
    gap = r_2.spos.values[0] - r_1.epos.values[0]
    if gap > max_gap:
        return False

    res = r_1.Ancestry.values[0] == r_2.Ancestry.values[0]
    return res


def merge_ancestries(r_1, r_2):
    """
    Replace the end position (physical and recombination)
    in r_1 with the values in r_2.
    Args:
        r_1: row 1
        r_2: row 2
    Returns: merged (1 row) pandas.DataFrame
    """
    merged = r_1.reset_index(drop=True)
    r_2 = r_2.reset_index(drop=True)
    merged.loc[0, 'epos'] = r_2.loc[0, 'epos']
    merged.loc[0, 'egpos'] = r_2.loc[0, 'egpos']

    return merged


def collapse_windows_to_tracts(df_chm):
    """
    If n continous windows are from the same ancestry,
    collpase to a single continous windows.
    Args:
        df_chm: pd.DataFrame
            Table with genomo coordinates and Ancestry assigments.
    Note:
        It is assumed that the input data represents one chromosome
        and only one Haplotype either maternal or paternal
    """
    # make sure df_chm is sorted by physical position
    df_chm = df_chm.sort_values(['spos'])

    current_ancstr_block = df_chm.iloc[[0]]
    tracts = list()

    for index in range(1, df_chm.shape[0]):
        current_row = df_chm.iloc[[index]]

        if is_ancestry_continuous(current_ancstr_block, current_row):
            current_ancstr_block = merge_ancestries(
                current_ancstr_block, current_row)

        else:
            tracts.append(current_ancstr_block)
            current_ancstr_block = current_row

    # add the last ancestry, the for loop does not include the last row

    tracts.append(current_ancstr_block)

    tracts = pd.concat(tracts).reset_index(drop=True)

    # compute tracts length
    tracts['len_bp'] = tracts.epos - tracts.spos
    tracts['len_cm'] = tracts.egpos - tracts.sgpos

    return tracts.reset_index(drop=True)


def _make_bins_to_fill_na(bins, ancestries):
    """
    Internal function to generate a table bins x ancestries
    """
    def myfmt(x, y): return str(round(x)) + '-' + str(round(y))
    bins_labs = [myfmt(a, b) for (a, b) in zip(bins[:-1], bins[1:])]
    df_1 = pd.DataFrame({'tract_length': bins_labs})
    df_2 = pd.DataFrame({'Ancestry': ancestries})
    return bins_labs, df_1.merge(df_2, how='cross')


def computer_tract_len_dist(bed_chr_h, nbins=50, max_cm_pos=250):
    """
    Computes the trackts length distribution for the given input.
    Note:
        It is assumed that the input data represents one chromosome
        and only one Haplotype either maternal or paternal
    Args:
        bed_chr_h: pd.DataFrame
            Table with data for only one chromosome and
            one Haplotype.
        nbins: int,
            the number of bins to use
        max_cm_pos: int,
            the maximum value for the binning)
    Returns:
        pd.DataFrame with distribution
    """
    # Collapse to continous ancestry tracts
    bed_chr_h = collapse_windows_to_tracts(bed_chr_h)
    bins = np.linspace(start=0, stop=max_cm_pos, num=nbins)

    # make useful label for the cm bins
    ancestries = bed_chr_h.Ancestry.unique().tolist()
    bins_labs, full_data = _make_bins_to_fill_na(bins, ancestries)

    bed_chr_h['tract_length'] = pd.cut(bed_chr_h.len_cm, bins,
                                       labels=bins_labs,
                                       include_lowest=True)

    distribution = (
        bed_chr_h[['Ancestry', 'tract_length']]
        .value_counts()
        .reset_index()
        .rename({0: 'Frequency'}, axis=1)
    )
    distribution['tract_length'] = distribution.tract_length.astype(str)
    distribution = pd.merge(full_data, distribution, how="left", on=[
                            'tract_length', 'Ancestry'])

    # Missing value indicates a frequency of 0
    return distribution.fillna(0)

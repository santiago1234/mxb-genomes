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
        fragments
    """

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


def collapse_windows_to_tracks(df_chm):
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
    tracks = list()

    for index in range(1, df_chm.shape[0]):
        current_row = df_chm.iloc[[index]]

        if is_ancestry_continuous(current_ancstr_block, current_row):
            current_ancstr_block = merge_ancestries(
                current_ancstr_block, current_row)

        else:
            tracks.append(current_ancstr_block)
            current_ancstr_block = current_row

    # add the last ancestry, the for loop does not include the last row

    tracks.append(current_ancstr_block)

    tracks = pd.concat(tracks).reset_index(drop=True)

    # compute tracks length
    tracks['len_bp'] = tracks.epos - tracks.spos
    tracks['len_cm'] = tracks.egpos - tracks.sgpos

    return tracks.reset_index(drop=True)


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
    bed_chr_h = collapse_windows_to_tracks(bed_chr_h)
    bins = np.linspace(start=0, stop=max_cm_pos, num=nbins)

    bed_chr_h['bins_cm'] = pd.cut(bed_chr_h.len_cm, bins,
                                  labels=bins.round()[:-1],
                                  include_lowest=True)
    new_names = {'index': 'tract_len_bin', 'bins_cm': 'relative_frequency'}
    distribution = (
        bed_chr_h
        .bins_cm
        .value_counts()
        .reset_index()
        .rename(new_names, axis=1)
    )

    return distribution

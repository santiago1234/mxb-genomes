# Plot tract length distribution

import numpy as np
import pandas as pd


def is_ancestry_continuous(r_1, r_2):
    """
    Params:
        r_1: pd.DataFrame, a single row pandas frame.
            this row contains the coulummn assigment
        r_2: same as r_1
    Returns:
        lgl: True if r_1 and r_2 represent the same ancestry
    Note:
        I assume that r_2 follows from r_1. I need
        to write conditions to make sure that r_1 and r_2
        if merged together represent continuous chromosome
        fragments
    """

    res = r_1.assignment.values[0] == r_2.assignment.values[0]
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
    """
    # make sure df_chm is sorted by physical position
    df_chm = df_chm.sort_values(['spos'])
    
    current_ancstr_block = df_chm.iloc[[0]]
    tracks = list()

    for index in range(1, df_chm.shape[0]):
        current_row = df_chm.iloc[[index]]

        if is_ancestry_continuous(current_ancstr_block, current_row):  
            current_ancstr_block = merge_ancestries(current_ancstr_block, current_row)

        else:
            tracks.append(current_ancstr_block)
            current_ancstr_block = current_row
            
    tracks = pd.concat(tracks).reset_index(drop=True)
    
    # compute tracks length
    tracks['len_bp'] = tracks.epos - tracks.spos
    tracks['len_cm'] = tracks.egpos - tracks.sgpos
    
    return tracks




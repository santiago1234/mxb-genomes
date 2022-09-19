"""
Script to compute fst in a set of input spectrums

usage:
    python scripts/fst_data.py {spectrum_files} {outfile.csv}
"""
import os
import sys
import pickle

import moments
import pandas as pd


sfile = 'data/joint-sfs-cat-missense-sim-10-pair-MXLxMXB.pkl'


def load_spectrum_file(sfile):
    """
    Load spectrum file
    """
    with open(sfile, 'rb') as data:
        spectrum = pickle.load(data)

    return spectrum


def fst(sfs) -> float:
    '''
    Compute fst
    '''
    assert len(sfs.shape) == 2, 'expected 2d spectrum'
    return sfs.fold().Fst()


def decode_file_name(sfile):
    '''
    decode the metadata in the input
    file name
    Return:
        tuple: (cat, sim-id, poppair)
    '''
    keywords = ['cat-', 'sim-', 'pair-']

    def get_key_word(key):
        """
        helper function to extract keyword
        """
        
        key_index = sfile.find(key)
        from_index_key = sfile[key_index:]
        keyword = from_index_key.split('-')[1]
        if keyword.endswith('.pkl'):
            return keyword[:-4]
        return keyword


    return [get_key_word(key) for key in keywords]


def main(sfiles: list, outtable: str):
    '''
    Calculates the fst for each spectrum in sfiles
    '''

    data = []
    for sfile in sfiles:
        metadata = decode_file_name(sfile)
        current_fst = fst(load_spectrum_file(sfile))
        metadata.append(current_fst)
        data.append(metadata)

    col_names = ['category', 'sim_id', 'pop_pair', 'fst']
    data = pd.DataFrame(data, columns=col_names)
    data.to_csv(outtable, index=False)



if __name__ == '__main__':

    *sfiles, outfile = sys.argv[1:]

    print(f'computing sfs for {len(sfiles)}')
    main(sfiles, outfile)

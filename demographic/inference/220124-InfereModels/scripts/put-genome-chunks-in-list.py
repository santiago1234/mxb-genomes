"""
Put together genome chunks in a python list object.
"""
import moments
import functools
import pickle
import glob
import gzip
from os import path
import re
import multiprocessing as mp
import sys

# Initialize global variables
PATH_TO_SPECTRUMS, PATH_TO_MLS, VARCAT, CORES, OUTFILE = sys.argv[1:]
CORES = int(CORES)
PROJECT_TO_SIZE = [30, 30, 30, 30]


def make_chunk_filenames(chunk_id):
    # Generate file names for a given category and chunk
    spec_file = f'spectrum_chunk_{chunk_id}_cat_{VARCAT}.pkl.gz'
    spec_file = path.join(PATH_TO_SPECTRUMS, spec_file)
    
    # mL files for intronic are named introns.
    vc = 'introns' if VARCAT == 'intronic' else VARCAT
    mL_file = f'mL_{vc}_chunk_{chunk_id}.txt'
    mL_file = path.join(PATH_TO_MLS, mL_file)

    return spec_file, mL_file


def read_spectrum_file(spec_file):
    '''Read moments.Spectrumb'''
    with gzip.open(spec_file, "rb") as f:
        sf = pickle.load(f)

    sf = sf.project(PROJECT_TO_SIZE)
    return sf


def read_mL_file(mL_file):
    '''Read mL from file'''
    f = open(mL_file, 'r')
    mL = f.readlines()[0]
    mL = mL.replace('mL:', '')
    mL = mL.strip()
    mL = float(mL)
    f.close()
    return mL


def load_chunk_data(chunk_id):
    sys.stderr.write(f'Loading data in chunk {chunk_id} ...\n')
    spec_file, mL_file = make_chunk_filenames(chunk_id)

    spec_chunk = read_spectrum_file(spec_file)
    mL_chunk = read_mL_file(mL_file)

    return (chunk_id, spec_chunk, mL_chunk)


if __name__ == '__main__':

    pattern = f'spectrum_chunk_*_cat_{VARCAT}.pkl.gz'
    fls = path.join(PATH_TO_SPECTRUMS, pattern)
    fls = glob.glob(fls)
    chunk_ids = [re.search('chunk_(\d*)_', x).group(1) for x in fls]

    pool = mp.Pool(CORES)
    results = pool.map(load_chunk_data, chunk_ids)
    pool.close()
    print('saving chunks')
    with gzip.open(OUTFILE, 'wb') as f:
        pickle.dump(results, f)
    print('done ...')

"""
usage:
    python construct-boostraped-replicates.py <PATH_TO_SPECTRUMS> <PATH_TO_MLS> <VARCAT> <OUTPUT_FILE>
Args:
    PATH_TO_SPECTRUMS: path to spectrum chunks
    PATH_TO_MLS: path to mLs chunks
    VARCAT: one of: synonymous, introns, intergenic, missense, or lof.
    OUTPUT_FILE: Filename for python object with boostrap replicates.

"""
import moments
import functools
import random
import pickle
import glob
import gzip
from os import path
import re
import sys


PATH_TO_SPECTRUMS, PATH_TO_MLS, VARCAT, OUT_REPLICATES = sys.argv[1:]

print(f'CATEGORY={VARCAT}')
# I assume the file name has the following format.
pattern = f'spectrum_chunk_*_cat_{VARCAT}.pkl.gz'
fls = path.join(PATH_TO_SPECTRUMS, pattern)
fls = glob.glob(fls)
chunk_ids = [re.search('chunk_(\d*)_', x).group(1) for x in fls]
print(f'processing {len(chunk_ids)} chunks.')


def make_chunk_filenames(chunk_id):
    # Generate file names for a given category and chunk
    spec_file = f'spectrum_chunk_{chunk_id}_cat_{VARCAT}.pkl.gz'
    spec_file = path.join(PATH_TO_SPECTRUMS, spec_file)
    
    # mL files for intronic are named introns.
    vc = 'introns' if VARCAT == 'intronic' else VARCAT
    mL_file = f'mL_{vc}_chunk_{chunk_id}.txt'
    mL_file = path.join(PATH_TO_MLS, mL_file)

    return spec_file, mL_file


def read_spectrum_file(spec_file):
    '''Read moments.Spectrumb'''
    with gzip.open(spec_file, "rb") as f:
        sf = pickle.load(f)

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

    print(f'Loading data in chunk {chunk_id} ...')
    spec_file, mL_file = make_chunk_filenames(chunk_id)

    spec_chunk = read_spectrum_file(spec_file)
    mL_chunk = read_mL_file(mL_file)

    return (chunk_id, spec_chunk, mL_chunk)


def get_genome_chunks(CHUNKS):
    '''
    Get the genome chunk data.
    '''
    CHUNKS = [load_chunk_data(i) for i in CHUNKS]
    return CHUNKS


def combine_spectrums(spec1, spec2):
    """
    Combine spectrum by adding the variants.
    Args:
        spec{1,2}: moments.Spectrum object
    """
    if len(spec1.shape) != len(spec2.shape):
        raise ValueError('Spectrums must have same dimension.')

    for x, y in zip(spec1.shape, spec2.shape):
        if x != y:
            print(x, y)
            raise ValueError(
                'Number of samples in each dimension should match.')

    for x, y in zip(spec1.pop_ids, spec2.pop_ids):
        if x != y:
            raise ValueError('Pops in each dimension must be the same.')

    return spec1 + spec2


def get_bootstrap_replicate(CHUNKS):
    """
    Args:
        CHUNKS: Genome chunks.
    Returns: (id, sfs, ml)
        id: string of the form 'j*k*l' where j, k, l are the chunk ids.
            It tell which data was used to generate the bootstraps.
        sfs: added spectrum
        ml: added ml
    """
    # sample with replacement
    boostrap_rep = random.choices(CHUNKS, k=len(CHUNKS))

    # Combine samples by adding them up
    boostrap_rep_sfs = functools.reduce(
        combine_spectrums, [x[1] for x in boostrap_rep])

    boostrap_rep_mL = sum([x[2] for x in boostrap_rep])

    boostrap_rep_id = '*'.join(str(x[0]) for x in boostrap_rep)

    print(f'Replicate: {boostrap_rep_id}')
    return boostrap_rep_id, boostrap_rep_sfs, boostrap_rep_mL


def get_boostraps(CHUNKS):
    '''There will be len(CHUNKS) boostraps replicates in output list.'''
    N = len(CHUNKS)

    boostraps = [get_bootstrap_replicate(CHUNKS) for _ in range(N)]
    return boostraps


#
# CODE
#
CHUNKS = get_genome_chunks(chunk_ids)
boostraps = get_boostraps(CHUNKS)
print('saving data')
boot_file = open(OUT_REPLICATES, 'wb')
pickle.dump(boostraps, boot_file)
boot_file.close()

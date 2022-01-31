"""
usage:
    python get-bootstrap-replicates.py CHUNKS.pkl.gz OUTFILE.pkl.gz CORES
"""
import sys
import multiprocessing
import functools
import random
import pickle
import gzip
import moments


CHUNKS_FILE, OUTFILE, CORES = sys.argv[1:]

CORES = int(CORES)


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
            sys.stderr.write(x, y)
            raise ValueError(
                'Number of samples in each dimension should match.')

    for x, y in zip(spec1.pop_ids, spec2.pop_ids):
        if x != y:
            raise ValueError('Pops in each dimension must be the same.')

    return spec1 + spec2


def get_bootstrap_replicate(CHUNKS, seed):
    """
    Args:
        seed: seed for random number generation
        CHUNKS: list
    Returns: (id, sfs, ml)
        id: string of the form 'j*k*l' where j, k, l are the chunk ids.
            It tell which data was used to generate the bootstraps.
        sfs: added spectrum
        ml: added ml
    """
    # sample with replacement
    random.seed(seed)
    boostrap_rep = random.choices(CHUNKS, k=len(CHUNKS))

    # Combine samples by adding them up
    boostrap_rep_sfs = functools.reduce(
        combine_spectrums, [x[1] for x in boostrap_rep])

    boostrap_rep_mL = sum([x[2] for x in boostrap_rep])

    boostrap_rep_id = '*'.join(str(x[0]) for x in boostrap_rep)

    sys.stderr.write(f'Replicate: {boostrap_rep_id}\n')
    return boostrap_rep_id, boostrap_rep_sfs, boostrap_rep_mL


if __name__ == '__main__':
    print('loading chunks ...')

    with gzip.open(CHUNKS_FILE, 'rb') as f:
        CHUNKS = pickle.load(f)
    print('Generating bootstrap replicates...')

    N_REPS = len(CHUNKS)
    pool = multiprocessing.Pool(CORES)

    BOOT_REP = pool.starmap(get_bootstrap_replicate, [
                            (CHUNKS, s) for s in range(N_REPS)])
    pool.close()
    print([x[0] for x in BOOT_REP])

    with gzip.open(OUTFILE, 'wb') as f:
        pickle.dump(BOOT_REP, f)
    print('done ...')

"""
Combined spectrum object into one.
Usage:
    python <output-file.pkl> <sp1> <sp2> <sp3> ... <spN>
"""
import sys
import functools
import pickle
import gzip
import moments

out_spectrum, *spec_files = sys.argv[1:]

# Functions
def load_spectrum(spec_gzip_file):
    """
    Load moments.Spectrum object from a pickled and
    gzziped file.
    """
    with gzip.open(spec_gzip_file, 'rb') as f:
        spectrum = pickle.load(f)
    return spectrum


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
            raise ValueError('Number of samples in each dimension should match.')

    for x, y in zip(spec1.pop_ids, spec2.pop_ids):
        if x != y:
            raise ValueError('Pops in each dimension must be the same.')

    return spec1 + spec2


# CODE
spectrums = [load_spectrum(x) for x in spec_files]
combined_spectrum = functools.reduce(combine_spectrums, spectrums)

outfile = open(out_spectrum, 'wb')
pickle.dump(combined_spectrum, outfile)
outfile.close()

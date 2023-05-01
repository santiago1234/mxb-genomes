"""
Get the 1d SFS from the joint SFS
"""
import pickle
import gzip

import moments
import pandas as pd

POPS =  ['CHB', 'IBS', 'MXB', 'YRI']

def load_ml(ml_file):
    """
    load ml
    """
    with open(ml_file, 'r') as file:
        number = float(file.read().strip())
    return number


def load_spectrum(spectrum_file):
    """
    load spectrum
    """
    with gzip.open(spectrum_file, "rb") as f:
        spectrum = pickle.load(f)

    return spectrum


def spectrum_fold_to_array(sf):
    '''Folds the expectrum and get a numpy array'''
    sf_folded = sf.fold()
    sf_folded = sf_folded[~sf_folded.mask].data
    return sf_folded


def to_frame(sfs_folded, pop):
    '''
    Put the folded SFS in a pandas frame, with metadata info
    '''
    minor_alle_f = list(range(1, len(sfs_folded)+1))

    d = {
        'Population': pop,
        'Frequency': sfs_folded,
        'Minor_allel_freq': minor_alle_f
    }
    return pd.DataFrame(d)


def single_sfs(popid, spectrum):
    """
    get single sfs
    """
    pops = spectrum.pop_ids.copy()

    indexes_to_marginalize = []

    for index, pop in enumerate(pops):
        if pop != popid:
            indexes_to_marginalize.append(index)

    sfs = spectrum.marginalize(indexes_to_marginalize)

    assert popid == sfs.pop_ids[0] 

    sfs = spectrum_fold_to_array(sfs)
    sfs = to_frame(sfs, popid)

    return sfs


def single_sfs_all_pops(spectrum):
    """
    get single sfs for all pops
    """
    sfs_list = []
    for pop in POPS:
        sfs = single_sfs(pop, spectrum)
        sfs_list.append(sfs)

    return pd.concat(sfs_list)


def marginalize(spectrum, poplist):
    """
    marginalize
    """
    pops = spectrum.pop_ids.copy()

    indexes_to_marginalize = []

    for index, pop in enumerate(pops):
        if pop not in poplist:
            indexes_to_marginalize.append(index)

    sfs = spectrum.marginalize(indexes_to_marginalize)

    return sfs


def Fst(spectrum):
    """
    Get the Fst from a joint spectrums
    """
    # get all pairwise SFS
    # first the pairwise poplist
    poppair = []
    for i in range(len(POPS)):
        for j in range(i+1, len(POPS)):
            poppair.append([POPS[i], POPS[j]])

    data = []
    # get the pairwise Fst

    for pair in poppair:
        sfs = marginalize(spectrum, pair)
        fst = sfs.Fst()
        data.append(pair + [fst])

    return pd.DataFrame(data, columns=['pop1', 'pop2', 'Fst'])


def quantile_data(quantile):
    """
    Loads the file in the quantile dir
    and gets the single SFS for each pop
    Also the Fst
    """
    spectrum_file = 'data/whole-genome/q{}/spectrum-cat_intergenic.pkl.gz'.format(quantile)
    ml_file = 'data/whole-genome/q{}/mL_intergenic.txt'.format(quantile)
    ml = load_ml(ml_file)
    spectrum = load_spectrum(spectrum_file)

    sfs = single_sfs_all_pops(spectrum)
    fst = Fst(spectrum)
    sfs['Quantile'] = quantile
    sfs['mL'] = ml
    fst['Quantile'] = quantile
    fst['mL'] = ml

    return sfs, fst


def main():
    quantiles = [1, 2, 3, 4]
    sfs_list = []
    sfss = pd.concat([quantile_data(q)[0] for q in quantiles])
    fsts = pd.concat([quantile_data(q)[1] for q in quantiles])
    sfss.to_csv('data/whole-genome/sfs-single-pop.csv', index=False)
    fsts.to_csv('data/whole-genome/fst.csv', index=False)


if __name__ == '__main__':
    main()

import pandas as pd
import glob
from os import path


def load_ind(ind_id, tracts_path, inter, end, haplo_names):
    """
    Returns a dict with the tracts data for two
    haplotypes.
    """
    names = ['chrn', 'spos', 'epos', 'Ancestry', 'sgpos', 'egpos']

    def load_file(hap):
        bedfile = ind_id + inter + hap + end + '.bed'
        bedfile = path.join(tracts_path, bedfile)
        return pd.read_csv(bedfile, names=names, sep='\t')

    hap_A = load_file(haplo_names[0])
    hap_A['haplo'] = 'A'
    hap_B = load_file(haplo_names[1])
    hap_B['haplo'] = 'B'

    return [hap_A, hap_B]


def load_tracts_pop(tracts_path, inter='_anc_', end='_cM', haplo_names=['A', 'B']):
    """
    Args:
        inter: string between individual label and haploid chromosome id in input file
        end: string at the end of input file. Note that the file should end in ".bed"
        haplo_name: The name for the haplotypes
    The file name for each tract file is formed by: ID + inter + haplo_names + end + '.bed'
    example: HG00554_anc_A_cM.bed

    Returns:
        a dict mapping individuals to a list [A, B] with the tracts data. Tracts data
        is represented as a pandas data frame.
    """
    bedfiles = glob.glob(tracts_path + '/' + '*.bed')
    ind_ids = [path.basename(x).split(inter)[0] for x in bedfiles]
    ind_ids = list(set(ind_ids))

    return {x: load_ind(x, tracts_path, inter, end, haplo_names) for x in ind_ids}

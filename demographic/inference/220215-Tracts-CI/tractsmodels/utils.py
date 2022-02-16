import os
import numpy as np

_defaul_labs = ['NAT', 'EUR', 'AFR']

LABELS = {
    'CLM': _defaul_labs,
    'MXL': _defaul_labs,
    'PEL': _defaul_labs,
    'PUR': ['AFR', 'EUR', 'NAT']
}


# string between individual label and haploid chromosome id in input file
inter = "_anc"
# string at the end of input file. Note that the file should end in ".bed"
end = "_cM.bed"


def list_individuals_in_dir(directory, inter='_anc', end='_cM.bed'):
    """
    Args:
        directory: Path to directory with tracts input bed files.
        inter: string between individual label and haploid chromosome id in input file
        end: string at the end of input file. Note that the file should end in ".bed"
    """
    _files = os.listdir(directory)
    files = [filex
             for filex in _files
             if filex.split('.')[-1] == "bed"]
    ind_names = list(set(file.split('_')[0] for file in files))
    return ind_names


def bootsamp(num):
    '''generates a list of positions of the samples to pick in a bootstrap'''
    return np.random.choice(range(num), replace=True, size=num)


def get_bootstrap_replicata(pop, seed):
    if seed == 0:
        return pop
    else:
        np.random.seed(seed)
        indivs = pop.indivs
        bootorder = bootsamp(len(indivs))
        indivs2 = [indivs[i] for i in bootorder]
        pop.indivs = indivs2
        return pop
        

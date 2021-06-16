import pandas as pd
import numpy as np
import sys
sys.path.append("../../")
from mxbgenomes.localancestry import tractsresults


# load the parameters t1 and t2

def load_tracts_result(subpop):
    """
    Load the tracts results for the given population
    Args:
        subpop: str one from MXL, PEL, CLM, PUR
    """
    dir_to_res = "output-" + subpop
    labels = ['NAT', 'EUR', 'AFR']
    # load the PUR data
    res = tractsresults.load_boot_rest(dir_to_res, boot_n=0, labels=labels)
    res['population'] = subpop
    # load the estimated parameters
    pars = tractsresults._load_boot_file(dir_to_res, filetype="pars")
    pars = pars.round(3)
    pars = "T1 = {}, T2 = {}".format(pars[0], pars[1])
    res['parameters'] = pars
    return res


for pop in ['PEL', 'CLM', 'PUR', 'MXL']:
    print(pop)
    tracts_res = load_tracts_result(pop)
    outname = "results/{}-res.csv".format(pop)
    tracts_res.to_csv(outname, index=False)

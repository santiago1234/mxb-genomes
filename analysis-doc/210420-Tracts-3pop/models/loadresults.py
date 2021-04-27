import pandas as pd
import os
import glob
import numpy as np
from scipy.stats import poisson


def _list_boostrap_files(dir_to_res, boot_n=0):
    """
    List the tracts output files for a particular boostrap run
    Returns:
        list of file for the given boostrap
    """
    start_pat = 'boot{n}_*'.format(n=boot_n)
    start_pat = os.path.join(dir_to_res, start_pat)
    return glob.glob(start_pat)


def _load_boot_file(dir_to_res, boot_n=0, filetype="dat"):
    """
    Args:
        filetype: str, either dat, bins, pred, etc.
    Loads the bins as an array
    """
    boot_file = [x for x in _list_boostrap_files(
        dir_to_res, boot_n) if x.endswith(filetype)][0]
    return np.loadtxt(boot_file)


def confIntMeanPoisson(mu, conf=0.8):
    return poisson.interval(alpha=conf, mu=mu + + 1e-9)


def load_boot_rest(dir_to_res, boot_n, labels):
    """
    Args:
        labels: list, the population labels, same order as
        it was used to fit the model
        boot_n: int, the boostrap number
    """
    bins = _load_boot_file(dir_to_res, boot_n, "bins")
    dat = _load_boot_file(dir_to_res, boot_n, "dat")
    pred = _load_boot_file(dir_to_res, boot_n, "pred")
    # Make data frame for dat

    def make_frame(arr_d, var_name):
        "helper function to construc data frame"
        df = pd.DataFrame(arr_d.T, columns=labels)
        df['bins'] = bins
        df = df.melt(id_vars=['bins'],
                     var_name="Ancestry", value_name=var_name)
        return df

    dat = make_frame(dat, 'dat')
    pred = make_frame(pred, 'pred')
    data = pd.merge(dat, pred, how='inner', on=['bins', 'Ancestry'])
    data['boot'] = boot_n
    c_i = data.pred.map(confIntMeanPoisson)
    data['ci_l'] = c_i.map(lambda x: x[0])
    data['ci_u'] = c_i.map(lambda x: x[1])
    return data



############
dir_to_res = "output/"
labels = ['NAT', 'EUR', 'AFR']
datos = load_boot_rest(dir_to_res, 0, labels)
datos.to_csv('mxl.csv', index=False)

import pandas as pd
import os
import glob
import numpy as np
import scipy.stats as st


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
    return data


def confIntMean(a, conf=0.95):
    mean, sem, m = np.mean(a), st.sem(a), st.t.ppf((1+conf)/2., len(a)-1)
    return mean, mean - m*sem, mean + m*sem

############
# TODO: put this into a function
dir_to_res = "output/"
labels = ['NAT', 'EUR', 'AFR']
datos = load_boot_rest(dir_to_res, 0, labels)

# compute mean and confidence interval

boots = list(range(1, 20))
boots = pd.concat([load_boot_rest(dir_to_res, x, labels) for x in boots])

# Compute the cofiden interval for preds
#
res = (
    boots.
    groupby(['bins', 'Ancestry'])['pred'].
    apply(confIntMean).
    reset_index()
)

res['mean'] = res.pred.map(lambda x: x[0])
res['ci_u'] = res.pred.map(lambda x: x[1])
res['ci_l'] = res.pred.map(lambda x: x[2])
res.drop(['pred'], axis=1, inplace=True)
datos = datos.merge(res, how='inner', on=['bins', 'Ancestry'])
datos.to_csv("mxl.csv", index=False)


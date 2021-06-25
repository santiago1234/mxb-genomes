import pandas as pd
import numpy as np
import sys
sys.path.append("../../")
from mxbgenomes.localancestry import tractsresults


def load_migmat(path_to_mdl_res):
    """
    load the migration matrix
    """
    migmat = [x for x in tractsresults._list_boostrap_files(
        path_to_mdl_res) if x.endswith('mig')]
    migmat = np.loadtxt(migmat[0])
    return migmat


def get_loglkl(path_to_mdl_res):
    """
    get the likelihood for the given model.
    """
    afile = [x for x in tractsresults._list_boostrap_files(
        path_to_mdl_res) if x.endswith('mig')][0]
    # a dirty way to get the number
    loglkl = afile[afile.find('boot0_') + 6: afile.find('_mig')]
    float(loglkl)
    return loglkl


def load_model_data(path_to_mdl_res, mdlname, labels):
    """
    Args:
        path_to_mdl_res: The path to the dir containing the tracts
            output files. It will load only the results for boostrap0.
        mdlname: the name of the model (e.g. ppx_xxp)
        labels: The population labels, same order at it was given in
            the fitting.
    """
    res = tractsresults.load_boot_rest(path_to_mdl_res,
            boot_n=0, labels=labels)
    res['model'] = mdlname
    res['loglkl'] = get_loglkl(path_to_mdl_res)

    # create the ancestry migration matrix.
    migmat = load_migmat(path_to_mdl_res)
    ancp = tractsresults.ancestry_prop_over_time(migmat, labels)
    ancp['model'] = mdlname

    # get the migration pulses.
    pulses = tractsresults.migration_pulses(migmat, labels)
    pulses['model'] = mdlname
    return res, ancp, pulses


ccxx_xxpp = load_model_data(
    "results/ccxx_xxpp/", "ccxx_xxpp", ['EUR', 'NAT', 'AFR', 'EAS'])
ccxx_xxpx_xxxp = load_model_data(
    "results/ccxx_xxpx_xxxp/", "ccxx_xxpx_xxxp", ['EUR', 'NAT', 'AFR', 'EAS'])
ppxx_ccxx_xxpx_xxxp = load_model_data("results/ppxx_ccxx_xxpx_xxpx/", 'ppxx_ccxx_xxpx_xxxp', ['EUR', 'NAT', 'AFR', 'EAS'])

fits = pd.concat([ccxx_xxpp[0], ccxx_xxpx_xxxp[0], ppxx_ccxx_xxpx_xxxp[0]])
ancp= pd.concat([ccxx_xxpp[1], ccxx_xxpx_xxxp[1], ppxx_ccxx_xxpx_xxxp[1]])
pulses= pd.concat([ccxx_xxpp[2], ccxx_xxpx_xxxp[2], ppxx_ccxx_xxpx_xxxp[2]])

fits.to_csv("results/fits.csv", index=False)
ancp.to_csv("results/ancp.csv", index=False)
pulses.to_csv("results/pulses.csv", index=False)

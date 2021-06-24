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


ppx_xxp = load_model_data("../210531-TRACTS-ppx_xxp/output-MXL/", "ppx_xxp", ['NAT', 'EUR', 'AFR'])
ppx_xxp_pxx = load_model_data("../210616-TRACTS-ppx_xxp_pxx/results", "ppx_xxp_pxx", ['EUR', 'NAT', 'AFR'])
ccx_xxp = load_model_data("../210619-TRACTS-ccx_xxp/results/", "ccx_xxp", ['EUR', 'NAT', 'AFR'])
ppx_ccx_xxp = load_model_data("../210619-TRACTS-ccx_xxp/ppx_ccx_xxp/results/", "ppx_ccx_xxp", ['EUR', 'NAT', 'AFR'])

fits = pd.concat([ppx_xxp[0], ppx_xxp_pxx[0], ccx_xxp[0], ppx_ccx_xxp[0]])
ancp = pd.concat([ppx_xxp[1], ppx_xxp_pxx[1], ccx_xxp[1], ppx_ccx_xxp[1]])
pulses = pd.concat([ppx_xxp[2], ppx_xxp_pxx[2], ccx_xxp[2], ppx_ccx_xxp[2]])

fits.to_csv("data/fits.csv", index=False)
ancp.to_csv("data/ancp.csv", index=False)
pulses.to_csv("data/pulses.csv", index=False)


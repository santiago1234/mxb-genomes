import pandas as pd
import numpy as np
import sys
sys.path.append("../../")
from mxbgenomes.localancestry import tractsresults


subpop = "MXL"
dir_to_res = "results/"
# same order as it was given to fit the data
labels = ['EUR', 'NAT', 'AFR']
# load the PUR data
res = tractsresults.load_boot_rest(dir_to_res, boot_n=0, labels=labels)
res['population'] = subpop

# load the estimated parameters
pars = tractsresults._load_boot_file(dir_to_res, filetype="pars")
pars = pars.round(4) * 100
pars = "init_Eu = {}, tstart = {}\nafam_time = {}, nuEu_time = {}".format(*pars)
res['parameters'] = pars

res.to_csv("results/results.csv", index=False)


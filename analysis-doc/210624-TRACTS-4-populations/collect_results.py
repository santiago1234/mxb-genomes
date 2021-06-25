import pandas as pd
import numpy as np
import sys
sys.path.append("../../")
from mxbgenomes.localancestry import tractsresults


dir_to_res = "results/ppxx_ccxx_xxpx_xxpx"
# same order as it was given to fit the data
labels = ['EUR', 'NAT', 'AFR', 'EAS']
# load the PUR data
res = tractsresults.load_boot_rest(dir_to_res, boot_n=0, labels=labels)

# load the estimated parameters
pars = tractsresults._load_boot_file(dir_to_res, filetype="pars")
pars = pars.round(4) * 100
print(pars)
res.to_csv("results/results.csv", index=False)


import sys

import pandas as pd

sys.path.append('../')
from mating import utils, randmating

tracts_path = '../../210514-ProcessXGMixOutForTracts/data/3-pops/tracts/MXL/'
poptracts = utils.load_tracts_pop(tracts_path)

newborn_id, hap_A, hap_B = randmating.newborn(poptracts)

outname = newborn_id + '_anc' + '_A_' + 'cM.bed'
# save the data
hap_A.to_csv(outname, header=False, index=False, sep='\t')
hap_B.to_csv(outname.replace('_A_', '_B_'), header=False, index=False, sep='\t')

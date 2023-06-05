import sys

import numpy as np
import pandas as pd

tracts_path = '../src/tracts-python3/'
tracs_models = '../../demographic/inference/220215-Tracts-CI/'

sys.path.append(tracts_path)
sys.path.append(tracs_models)

import tracts
from tractsmodels import models, utils
# note the namming is wrong
import ppxx_ccxx_xxpx_xxxp


BOOTNUM = sys.argv[1]  # Pass an integer value
BOOTNUM = int(BOOTNUM)
POPULATION = 'PEL'
MODEL = 'ppxx_ccxx_xxpx_xxxp'
dir_to_data = '../210514-ProcessXGMixOutForTracts/data/4-pops/tracts/PEL/'

startparams_ppxx_ccxx_xxpx_xxxp = [
    0.9,   # initEu_frac
    0.01,  # frac1
    0.12,  # frac2
    0.145,  # t1
    0.054,  # frac3
    0.014,  # frac4
    0.08,  # t2
    0.05,  # t3
]

MODELS = {
    'ppxx_ccxx_xxpx_xxxp': (ppxx_ccxx_xxpx_xxxp.ppxx_ccxx_xxpx_xxxp, ppxx_ccxx_xxpx_xxxp.outofbounds_ppxx_ccxx_xxpx_xxxp, startparams_ppxx_ccxx_xxpx_xxxp)
}


labels = ['EUR', 'NAT', 'AFR', 'EAS']

# ----------------------------------------------------------------

func, bound, startparams = MODELS[MODEL]
outf = f'results/inference-PEL-4pops/{MODEL}-boot{BOOTNUM}'

inter = "_anc"
end = "_cM.bed"

names = utils.list_individuals_in_dir(dir_to_data, inter, end)

pop = tracts.population(names=names, fname=(dir_to_data, inter, end))

# generate a boostrap replicate
# when BOOTNUM is 0 we get back the same data (not a resample)
pop = utils.get_bootstrap_replicata(pop, BOOTNUM)

(bins, data) = pop.get_global_tractlengths(npts=50)
data = [data[poplab] for poplab in labels]
# Calculate ancestry proportions
props = list(pop.get_mean_ancestry_proportions(labels))
Ls = pop.Ls
nind = pop.nind
cutoff = 2  # ignore particular bins


xopt = tracts.optimize_cob_fracs2(startparams, bins, Ls, data, nind,
                                  func, props, outofbounds_fun=bound, cutoff=cutoff, epsilon=1e-2)


print(xopt)

optmod = tracts.demographic_model(func(xopt, props))

optpars = xopt
maxlik = optmod.loglik(bins, Ls, data, pop.nind, cutoff=cutoff)

expects = []
for popnum in range(len(data)):
    expects.append(optmod.expectperbin(Ls, popnum, bins))

expects = nind * np.array(expects)

bootnum = 0

# save the poplabels order and the likelihood
info = pd.DataFrame(list(enumerate(labels)), columns=['order', 'Anc'])
info['loglik'] = maxlik
info['mdl'] = MODEL
info['Population'] = POPULATION
info['BOOT'] = BOOTNUM
info.to_csv(outf + '_info.tsv', sep='\t', index=False)


fbins = open(outf + "_bins", 'w')
fbins.write("\t".join(map(str, bins)))
fbins.close()


fdat = open(outf + "_dat", 'w')
for popnum in range(len(data)):
    fdat.write("\t".join(map(str, data[popnum])) + "\n")

fdat.close()
fmig = open(outf + "_mig", 'w')
for line in optmod.mig:
    fmig.write("\t".join(map(str, line)) + "\n")

fmig.close()
fpred = open(outf + "_pred", 'w')
for popnum in range(len(data)):
    fpred.write(
        "\t".join(map(str, pop.nind * np.array(optmod.expectperbin(Ls, popnum, bins)))) + "\n")
fpred.close()


fpars = open(outf + "_pars", 'w')
fpars.write("\t".join(map(str, optpars)) + "\n")
fpars.close()

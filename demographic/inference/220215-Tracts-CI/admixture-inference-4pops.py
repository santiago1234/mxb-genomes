import numpy as np
import pandas as pd
from tractsmodels import models, utils
import sys
tracts_path = '../../../analysis-doc/src/tracts-python3'
sys.path.append(tracts_path)
path_to_4pops_mdls = '../../../analysis-doc/210624-TRACTS-4-populations/models'
sys.path.append(path_to_4pops_mdls)
import tracts
import ccxx_xxpp, ppxx_ccxx_xxpp

MODEL, BOOTNUM = sys.argv[1:]
BOOTNUM = int(BOOTNUM)
dir_to_data = '../../../analysis-doc/210514-ProcessXGMixOutForTracts/data/4-pops/tracts/MXL/'
POPULATION = 'MXL'
# ----------------------------------------------------------------
# start parameters and models

startparams_ccxx_xxpp = [
    0.05,  # frac1
    0.05,  # frac2
    0.158,  # t1
    0.067,  # frac3
    0.007,  # frac4
    0.116  # t2
]

startparams_ppxx_ccxx_xxpp = [
    0.238,   # initEu_frac
    0.058,  # frac1
    0.049,  # frac2
    0.159,  # t1
    0.113,  # frac3
    0.007,  # frac4
    0.128,  # t2
]


MODELS = {
    'ccxx_xxpp': (ccxx_xxpp.ccxx_xxpp, ccxx_xxpp.outofbounds_ccxx_xxpp, startparams_ccxx_xxpp),
    'ppxx_ccxx_xxpp': (ppxx_ccxx_xxpp.ppxx_ccxx_xxpp, ppxx_ccxx_xxpp.outofbounds_ppxx_ccxx_xxpp, startparams_ppxx_ccxx_xxpp)
}

labels = ['EUR', 'NAT', 'AFR', 'EAS']
# ----------------------------------------------------------------

func, bound, startparams = MODELS[MODEL]
outf = f'results/inference-MXL-4pops/{MODEL}-boot{BOOTNUM}'


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

if MODEL in ['ccx_xxp', 'ppx_ccx_xxp']:
    # I do this because we fix the ancestry
    # proportions in the ppx_xxp and ppx_xxp_pxx models
    optmod = tracts.demographic_model(func(xopt))
else:
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

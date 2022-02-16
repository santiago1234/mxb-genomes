import numpy as np
import pandas as pd
from tractsmodels import models, utils
import sys
tracts_path = '../../../analysis-doc/src/tracts-python3'
sys.path.append(tracts_path)
import tracts

pops = ['PEL', 'MXL', 'CLM', 'PUR']
dir_to_data = '../../../analysis-doc/210514-ProcessXGMixOutForTracts/data/3-pops/tracts/PEL/'
MODEL = 'ppx_xxp'
BOOTNUM = 0

# Infere population name from dir_to_data
POPULATION, *_ = [x for x in pops if x in dir_to_data] 
func, bound, startparams = models.MODELS[MODEL]
labels = utils.LABELS['PEL']

outf = f'results/inference/{POPULATION}-{MODEL}-boot{BOOTNUM}'

inter = "_anc"
end = "_cM.bed"

names = utils.list_individuals_in_dir(dir_to_data, inter, end)

pop = tracts.population(names=names, fname=(dir_to_data, inter, end))

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

if MODEL == 'ppx_xxp':
    # I do this because we fix the ancestry
    # proportions in the ppx_xxp model
    optmod = tracts.demographic_model(func(xopt, props))
else:
    optmod = tracts.demographic_model(func(xopt))


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

#!/usr/bin/env python

import os
import sys
import scipy
tractspath = "../../src/tracts-python3"  # the path to tracts if not in your default pythonpath
sys.path.append(tractspath)
import tracts
import ppx_xxp_pxx_fix_tstart
import numpy
import pylab

from warnings import warn

# optimization method. Here we use brute force; other options are implemented
# in tracts but are not implemented in this driver script.
method = "brute"

# demographic models to use

# ppx_xxp_fix has an initial pulse of migration from populations 0 and 1,
# followed by a pulse from population 2. "fix" referes to the fact that the
# migration rates in the model are fixed to the observed global ancestry
# proportions--then, we only have to optimize the timing of the migrations.
func = ppx_xxp_pxx_fix_tstart.ppx_xxp_pxx_fix_tstart

# this function keeps track of whether parameters are in a "forbidden" region:
# whether mproportions are between 0 and 1, times positive, etc.
bound = ppx_xxp_pxx_fix_tstart.outofbounds_ppx_xxp_pxx_fix_start

# defines the values of parameters to loop over in the brute force
# step. Times are in units of 100 generations.
slices = (
        slice(.1, .5, .02),  # range to explore for init_Eu
        slice(0.06, .1, .01),  # " afam_time
        slice(0.06, .9, .02)   # " nuEu_time
    )


# absolute bounds that parameters are not allowed to cross: times
# must be between 0 and 100 generations.
bounds = [(0, 1), (0, 1)]

# choose order of populations: the labels in our ancestry files are
# strings, and we need to tell tracts which string correspond to which
# population in the model. Here the population labels are (somewhat
# confusingly) numbers that do not match the order in the population.
# "labels" will tell tracts that model population 0 has label '0', model
# population 1 has label '2', and model population 2 has label '1' in the
# local ancestry files.
labels = ['EUR', 'NAT', 'AFR']

# directories in which to read and write
directory = "../../210514-ProcessXGMixOutForTracts/data/3-pops/tracts/MXL/"
outdir = "results/"

if not os.path.exists(outdir):
    os.makedirs(outdir)


# string between individual label and haploid chromosome id in input file
inter = "_anc"
# string at the end of input file. Note that the file should end in ".bed"
end = "_cM.bed"


# Get a list of all individuals in directory.
_files = os.listdir(directory)
files = [file
        for file in _files
        if file.split('.')[-1] == "bed"]  # only consider bed files

# Get unique individual labels
names = list(set(file.split('_')[0] for file in files))

if len(_files) != len(files):
    warn("some files in the bed directory were ignored, since they do not "
            "end with `.bed`.")

# Load the population using the population class's constructor. It
# automatically iterates over individuals and haploid copies (labeled _A"
# and "_B" by default
pop = tracts.population(names=names, fname=(directory, inter, end))

# Rather than creating a new population for each bootstrap instance, we
# just replace the list of individuals to iterate over. We need to save a
# copy of the initial list of individuals to do this!
indivs = pop.indivs

# generate the histogram of tract lengths
(bins, data) = pop.get_global_tractlengths(npts=50)

data = [data[poplab] for poplab in labels]

# Calculate ancestry proportions
props = list(pop.get_mean_ancestry_proportions(labels))

Ls = pop.Ls
nind = pop.nind

cutoff = 2  # ignore particular bins

bootorder = range(len(indivs))


xopt = tracts.optimize_brute_fracs2(
    bins, Ls, data, nind, func, props, slices, outofbounds_fun=bound, cutoff=cutoff)
print(xopt)
optmod = tracts.demographic_model(func(xopt[0], props))
optpars = xopt[0]
liks = xopt[1]
maxlik = optmod.loglik(bins, Ls, data, pop.nind, cutoff=cutoff)

expects = []
for popnum in range(len(data)):
    expects.append(optmod.expectperbin(Ls, popnum, bins))

expects = nind * numpy.array(expects)

bootnum = 0
outf = outdir + "boot%d_%2.2f" % (bootnum, maxlik,)

fbins = open(outf + "_bins", 'w')
fbins.write("\t".join(map(str, bins)))
fbins.close()

fliks = open(outf + "_liks", 'w')
fliks.write("\t".join(map(str, liks)))
fliks.close()

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
        "\t".join(map(str, pop.nind * numpy.array(optmod.expectperbin(Ls, popnum, bins)))) + "\n")
fpred.close()

fpars = open(outf + "_pars", 'w')

fpars.write("\t".join(map(str, optpars)) + "\n")
fpars.close()

# bootstrap order in case we need to rerun/check something.
ford = open(outf + "_ord", 'w')
ford.write("\t".join(map(lambda i: "%d" % (i,), bootorder)))
ford.close()

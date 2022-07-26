import sys
import fwdpy11
import time
import gzip
sys.path.append('../')
from simutils import utils

# Load the test data

sim_dat = utils.simuldata(path_to_samples='test-data/',
                          sample_id=23, path_to_genetic_maps='test-data/')
random_seed = sys.argv[1:][0]
random_seed = int(random_seed)
out_file = f'results/simulations/sim-seed-{random_seed}-pop.bin'
print(f'simulation for: {out_file}')


mut_labels = {
    'neutral': 0,
    'missense': 1,
    'synonymous': 2,
    'LOF': 3,
}

# Make neutral geneic regions


nregions = []
for _, noexon in sim_dat.noncoding_intervals.iterrows():
    nregions.append(
        fwdpy11.Region(beg=noexon.start, end=noexon.end,
                       weight=sim_dat.m_noncoding, label=mut_labels['neutral'])
    )

# synonymous we assume they are neutral
for _, exon in sim_dat.coding_intervals.iterrows():
    nregions.append(
        fwdpy11.Region(beg=exon.start, end=exon.end,
                       weight=sim_dat.m_synonymous, label=mut_labels['synonymous'])
    )

# DFEs and selected regions
Ne = 11372.91
shape = 0.1596
scale = 2332.3
mean_s = (shape * scale) / (2 * Ne)


# DFE for LOF
shape_lof = 0.3589
scale_lof = 7830.5
mean_s_lof = (shape_lof * scale_lof) / (2 * Ne)


sregions = []
for _, exon in sim_dat.coding_intervals.iterrows():
    # missense
    sregions.append(
        fwdpy11.GammaS(
            beg=exon.start, end=exon.end,
            weight=sim_dat.m_missense,
            mean=-mean_s,  # NOTE: Mean is negative
            shape_parameter=shape,
            h=1,
            label=mut_labels['missense'])
    )
    # loss of function
    sregions.append(
        fwdpy11.GammaS(
            beg=exon.start, end=exon.end,
            weight=sim_dat.m_LOF,
            mean=-mean_s_lof,  # NOTE: Mean is negative
            shape_parameter=shape_lof,
            h=1,
            label=mut_labels['LOF'])
    )

# Recombination regions
nrec = len(sim_dat.rmap) - 1
recregions = []
for i in range(nrec):
    recregions.append(
        fwdpy11.PoissonInterval(
            beg=sim_dat.rmap.left[i],
            end=sim_dat.rmap.right[i],
            mean=sim_dat.rmap.rate[i] * sim_dat.rmap.span[i]
        )
    )

# Rates

# The neutral mutation rate and selected mutation rate
neutral_ml = sim_dat.ml_noncoding + sim_dat.ml_synonymous
selected_ml = sim_dat.ml_missense + sim_dat.ml_LOF

rates = fwdpy11.MutationAndRecombinationRates(
    neutral_mutation_rate=neutral_ml,
    selected_mutation_rate=selected_ml,
    recombination_rate=None)


Ne = 10000
pop = fwdpy11.DiploidPopulation(N=Ne, length=int(1e6))
pop.N
pop.tables.genome_length


SIM_LEN = 10 * pop.N

# the parameters that fwdpy11 needs to run the simulation
p = {
    # neutral mutations (none for now, can add after the fact)
    "nregions": nregions,
    "gvalue": fwdpy11.Multiplicative(2.0),  # fitness model
    "sregions": sregions,
    "recregions": recregions,
    "rates": rates,
    "prune_selected": True,
    "demography": fwdpy11.DiscreteDemography(),  # pass the demographic model
    "simlen": SIM_LEN
}
params = fwdpy11.ModelParams(**p)
# run the simulation
# set up the random number generator
rng = fwdpy11.GSLrng(random_seed)

# run the simulation
print('runnning simulation ...')
time1 = time.time()
fwdpy11.evolvets(
    rng, pop, params, simplification_interval=100, suppress_table_indexing=True
)
print("Simulation took", int(time.time() - time1), "seconds")

# simulation finished
print("Final population sizes =", pop.deme_sizes())

pop.dump_to_file(out_file)

import sys
import fwdpy11
import demes
import time
import gzip
sys.path.append('../')
from simutils import utils
from simutils.utils import DFE_missense, DFE_lof


## Demographic model

graph = '../220422-fwdpy11-initial-test/ADMIXTURE-MXL.yml'
graph = demes.load(graph)
demog = fwdpy11.discrete_demography.from_demes(graph)

print('demography set ..')

# Load the test data

sim_dat = utils.simuldata(path_to_samples='../220506-SetSimulationDesign/test-data/',
                          sample_id=23, path_to_genetic_maps='../220506-SetSimulationDesign/test-data/')

random_seed = int(42)
out_file = f'test-sim-pop.bin'
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

## SELECTED REGIONS

sregions = []
for _, exon in sim_dat.coding_intervals.iterrows():
    # missense
    sregions.append(
        fwdpy11.GammaS(
            beg=exon.start, end=exon.end,
            weight=sim_dat.m_missense,
            mean=-DFE_missense.mean(),  # NOTE: Mean is negative
            shape_parameter=DFE_missense.shape,
            h=1,
            label=mut_labels['missense'])
    )
    # loss of function
    sregions.append(
        fwdpy11.GammaS(
            beg=exon.start, end=exon.end,
            weight=sim_dat.m_LOF,
            mean=-DFE_missense.mean(),  # NOTE: Mean is negative
            shape_parameter=DFE_missense.shape,
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


## set the population

print("Initial sizes =", demog.metadata["initial_sizes"])
initial_sizes = [
    demog.metadata["initial_sizes"][i]
    for i in sorted(demog.metadata["initial_sizes"].keys())
]
pop = fwdpy11.DiploidPopulation(initial_sizes, sim_dat.end - sim_dat.start)


print(f'Pop info N={pop.N}, genome length = {pop.tables.genome_length}')


SIM_LEN = 10 * pop.N
print(f'Simulation will run for {SIM_LEN} generations.')

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

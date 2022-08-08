import sys
import fwdpy11
import demes
import time
import gzip
sys.path.append('../')
from simutils import utils
from simutils.utils import DFE_missense, DFE_lof


sim_id = sys.argv[1:][0]
## Demographic model

graph = '../220422-fwdpy11-initial-test/ADMIXTURE-MXL.yml'
graph = demes.load(graph)
demog = fwdpy11.discrete_demography.from_demes(graph)

print('demography set ..', file = sys.stderr)

# Load the test data

sim_dat = utils.simuldata(path_to_samples='../../data/220404-SimulationData/data/samples/',
                          sample_id=sim_id, path_to_genetic_maps='../../../resources/genetic-maps/')

print('Starting simulation for:', file = sys.stderr)
print(sim_dat)

random_seed = int(42)
out_file = f'results/simulations/sim-{sim_id}-pop.bin'
print(f'output file is: {out_file}', file = sys.stderr)


mut_labels = {
    'neutral': 0,
    'missense': 1,
    'synonymous': 2,
    'LOF': 3,
}


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

print("Initial sizes =", demog.metadata["initial_sizes"], file = sys.stderr)
initial_sizes = [
    demog.metadata["initial_sizes"][i]
    for i in sorted(demog.metadata["initial_sizes"].keys())
]

# genome length
L = sim_dat.end - sim_dat.start

pop = fwdpy11.DiploidPopulation(initial_sizes, L)


print(f'Pop info N={pop.N}, genome length = {pop.tables.genome_length}', file = sys.stderr)


# the parameters that fwdpy11 needs to run the simulation
p = {
    # neutral mutations (none for now, can add after the fact)
    "nregions": [],
    "gvalue": fwdpy11.Multiplicative(2.0),  # fitness model
    "sregions": sregions,
    "recregions": recregions,
    "rates": rates,
    "prune_selected": True,
    "demography": demog,
    "simlen": demog.metadata["total_simulation_length"],  # the total time to simulate
}

params = fwdpy11.ModelParams(**p)
# run the simulation
# set up the random number generator
rng = fwdpy11.GSLrng(random_seed)

# run the simulation
print('runnning simulation ...', file = sys.stderr)
print(f'Simulation lenght: {demog.metadata["total_simulation_length"]}')

time1 = time.time()
fwdpy11.evolvets(
    rng, pop, params, simplification_interval=100, suppress_table_indexing=True
)
print("Simulation took", int(time.time() - time1), "seconds", file = sys.stderr)

# simulation finished
print("Final population sizes =", pop.deme_sizes(), file = sys.stderr)

pop.dump_to_file(out_file)

import fwdpy11
import demes
import gzip
import time


# load the graph and import it as a fwdpy11 demography object
graph = demes.load('ADMIXTURE-MXL.yml')
demog = fwdpy11.discrete_demography.from_demes(graph)

# set up the region to simulate
L = 10 * 1e6  # length in base pairs
r = 1e-8  # the recombination rate per base pair
print("r * L =", r * L)


# initialize the ancestral population
print("Initial sizes =", demog.metadata["initial_sizes"])
initial_sizes = [
    demog.metadata["initial_sizes"][i]
    for i in sorted(demog.metadata["initial_sizes"].keys())
]
pop = fwdpy11.DiploidPopulation(initial_sizes, L)


# set up the random number generator
rng = fwdpy11.GSLrng(54321)  # should take a seed in the future

# the parameters that fwdpy11 needs to run the simulation
p = {
    "nregions": [],  # neutral mutations (none for now, can add after the fact)
    "gvalue": fwdpy11.Additive(2.0),  # fitness model
    "sregions": [],  # selected mutations (none for now)
    "recregions": [fwdpy11.PoissonInterval(0, L, r * L)],  # total rec rate in region
    "rates": (0.0, 0.0, None),  # neutral and selected mutations have rates 0
    "prune_selected": False,
    "demography": demog,  # pass the demographic model
    "simlen": demog.metadata["total_simulation_length"],  # the total time to simulate
}
params = fwdpy11.ModelParams(**p)


# run the simulation
print('runnning simulation ...')
time1 = time.time()
fwdpy11.evolvets(
    rng, pop, params, simplification_interval=100, suppress_table_indexing=True
)
print("Simulation took", int(time.time() - time1), "seconds")

# simulation finished
print("Final population sizes =", pop.deme_sizes())

# add neutral mutations
# The mutation rate, per haploid genome per generation
print('adding neutral mutations')
u = 1e-8
nmuts = fwdpy11.infinite_sites(rng, pop, u*L)
print(f"{nmuts} neutral mutations added")

# save the simulation results
with gzip.open('data/sim-pop.gz', 'wb') as f:
    pop.pickle_to_file(f)

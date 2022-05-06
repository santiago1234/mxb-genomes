import msprime
import demes
import tskit
import gzip

# Simulation parameters. This parameters
# should match those in sim-fwdpy11.py
rec_rate = 1e-8
# the mutation rate
u = 1e-8
N_ind_per_deme = 25  #samples per population 
L = 10 * 1e6


graph = demes.load("ADMIXTURE-MXL.yml")

print('setting demography')
demography = msprime.Demography.from_demes(graph)

print('setting evolution model')
evolution_model = [msprime.DiscreteTimeWrightFisher(duration=20), msprime.StandardCoalescent()]

# Samples to simulate
demes = ['YRI', 'IBS', 'CHB', 'MXB', 'MXL']
Samples = [msprime.SampleSet(N_ind_per_deme, population=x) for x in demes]

print('running simulation ...')
ts = msprime.sim_ancestry(
    samples=Samples,
    demography=demography,
    recombination_rate=rec_rate, 
    sequence_length=L,
    model=evolution_model,
    ploidy=2,
    random_seed=42
)

print('adding mutations ...')
# We set discrete_genome=False to use the infinite sites model
mts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)


mts.dump('data/ts-msprime.ts')

with open("data/simulated-genomes.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file)

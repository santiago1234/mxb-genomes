import pandas as pd
import msprime
import demes
import demesdraw
import tskit


print('loaded packages')

graph = demes.load("Best-Model-Intronic.yml")
demography = msprime.Demography.from_demes(graph)

print('setting demography')

evolution_model = [msprime.DiscreteTimeWrightFisher(duration=20), msprime.StandardCoalescent()]

print('setting evolution model')

## Rate Map
chr22map = pd.read_csv(
    '/data/users/smedina/data-resources/genetics-maps/beagle-genetic_maps/GRCh38-genetic-map/plink.chr22.GRCh38.map',
    sep=' ', names = ['chr', '?', 'rate', 'position']
)


position = chr22map['position'].tolist()
position.insert(0, 0)
rate = chr22map['rate'].tolist()
rate_map = msprime.RateMap(position=position, rate=rate)

print('setting genetic map')

## Simulation
N = 20
Samples = [
    msprime.SampleSet(N, population='YRI'),
    msprime.SampleSet(N, population='IBS'),
    msprime.SampleSet(N, population='CHB'),
    msprime.SampleSet(N, population='MXB')
]

print('running simulation ...')

ts = msprime.sim_ancestry(
    samples=Samples,
    demography=demography,
    recombination_rate=rate_map, 
    model= evolution_model,
    ploidy=2,
    random_seed=42
)

print('adding mutations ...')

mts = msprime.sim_mutations(ts, rate=1e-6, random_seed=42)


print('saving data ...')

pops = ['YRI']*20 + ['IBS']*20 + ['CHB']*20 + ['MXB']*20
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"{pops[i]}_{str(i)}indv" for i in range(n_dip_indv)]

with open("simulated-genomes-chr22.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file, individual_names=indv_names)


import pandas as pd
import msprime
import demes
import demesdraw
import tskit


print('loaded packages')

graph = demes.load("../../inference/220225-Combine-Inferences/ADMIXTURE.yml")
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
N_YRI = 108
N_IBS = 107
N_CHB = 103
N_MXB = 50
N_MXL = 64
N_PEL = 85
N_CLM = 94
N_PUR = 104

Samples = [
    msprime.SampleSet(N_YRI, population='YRI'),
    msprime.SampleSet(N_IBS, population='IBS'),
    msprime.SampleSet(N_CHB, population='CHB'),
    msprime.SampleSet(N_MXB, population='MXB'),
    msprime.SampleSet(N_MXL, population='MXL'),
    msprime.SampleSet(N_PEL, population='PEL'),
    msprime.SampleSet(N_CLM, population='CLM'),
    msprime.SampleSet(N_PUR, population='PUR')
]

print('running simulation ...')

ts = msprime.sim_ancestry(
    samples=Samples,
    demography=demography,
    recombination_rate=1e-8, 
    sequence_length=50000000, # 100Mb
    model= evolution_model,
    ploidy=2,
    random_seed=42
)

print('adding mutations ...')

mts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)

print('saving data ...')

pops = ['YRI']*N_YRI + ['IBS']*N_IBS + ['CHB']*N_CHB + ['MXB']*N_MXB + ['MXL']*N_MXL + ['PEL']*N_PEL + ['CLM']*N_CLM + ['PUR']*N_PUR
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"{pops[i]}_{str(i)}indv" for i in range(n_dip_indv)]

with open("data/simulated-genomes-chr22.vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file, individual_names=indv_names)


import moments
import demes
import pandas as pd

print('loading model')
graph = demes.load('ADMIXTURE-MXL.yml')

# I use this dic just to have the same names as in the other simulations
demes_ = {
    4: 'IBS',
    6: 'MXB',
    7: 'MXL'
}
inv_map = {v: k for k, v in demes_.items()}


def sfs_single(pop_id, N):
    """
    Args:
        apop: DiploidPopulation
        pop_id: the deme id
        N: number of nodes (diploid number) to include in the computation
    """
    sampled_demes=[pop_id]
    sf = moments.Spectrum.from_demes(
        graph,
        sampled_demes=sampled_demes,
        sample_sizes=[2*N]
    )
    print(f'pod id: {pop_id}')
    sf = sf.data
    return pd.DataFrame(
        {'F': sf,
         'derived_allel_freq': range(len(sf)),
         'pop_id': inv_map[pop_id]
        }
    )


print('getting SFS')
sfs = [sfs_single(demes_[i], 25) for i in demes_.keys()]
pd.concat(sfs).assign(Simulator='moments').to_csv('results/simulated-sfss-moments.csv', index=False)

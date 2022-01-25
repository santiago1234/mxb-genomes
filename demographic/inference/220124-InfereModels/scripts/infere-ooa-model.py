'''
usage:
    python scripts/infere-ooa-model.py <model> <parameters> <mL_file> <sfs_data> <outprefix>

Args:
    <model>: A YAML file in ``demes`` format.
    <parameters> inference_options:
    <mL_file> file with the mutation rate
    <sfs_data> The SFS to fit, which must have pop_ids specified.
    <outprefix> File prefix for output files.
'''
import sys
import pandas as pd
import moments
import pickle
import gzip


model, parameters, mL_file, sfs_data, outprefix = sys.argv[1:]
best_model = outprefix + '.yml'
best_parameters = outprefix + '-bestparameters.csv'

PROJECT_TO_SIZE = [40, 40, 40]

def read_mL_from_file(mL_file):
    '''Read mL from file'''
    f = open(mL_file, 'r')
    mL = f.readlines()[0]
    mL = mL.replace('mL:', '')
    mL = mL.strip()
    mL = float(mL)
    f.close()
    return mL


mL = read_mL_from_file(mL_file)
print(f'mL: {mL}')


print('Loading data...')
with gzip.open(sfs_data, "rb") as f:
    sf = pickle.load(f)

print(sf.pop_ids)

# MARGINALIZE AND FOLD
sf = sf.marginalize([sf.pop_ids.index('MXB')])

print('Projecting and folding data...')
sf = sf.project(PROJECT_TO_SIZE)
sf = sf.fold()


ret = moments.Demes.Inference.optimize(
    deme_graph=model,
    output=best_model,
    inference_options=parameters,
    data=sf,
    verbose=50,
    fit_ancestral_misid=False,
    uL=mL,
    overwrite=True
)


# Save best parameters in a table
pd.DataFrame({'parameter': ret[0], 'value': ret[1]}).to_csv(best_parameters, index=False)

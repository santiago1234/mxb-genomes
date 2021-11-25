"""
usage:
    python inference.py <model.yml> <options.yml> <fitted.yml> <outparams.txt>
"""
import sys
import moments
import demes
import demesdraw
import matplotlib.pyplot as plt
import pickle
import gzip

model_yml, options_params, inferred_model_yml, outfile_lkl = sys.argv[1:]
model_yml = "ooa-4pops.yml"
options_params = "ooa4pos_options.yaml"
inferred_model_yml = "results/ooa-4pops-fitted.yml"

# data
sf = "../../data/210804-Compute-jSFS/data/spectrums/5d-csq-synonymous-spectrum.pkl.gz"
with gzip.open(sf, "rb") as f:
    sf = pickle.load(f)

print(sf.pop_ids)

## marginalize
sf = sf.marginalize([3])
print("marginalizing to: ", sf.pop_ids)


## project and fold
n = 35
sf = sf.project([n] * 4)
## 
sf = sf.fold()
print(sf.shape)

## fitting
print('starting optimization ...')
uL = 0.14419746897690008
ret = moments.Demes.Inference.optimize(
    deme_graph=model_yml,
    output=inferred_model_yml,
    inference_options=options_params,
    data=sf,
    verbose=10,
    maxiter=1000,
    uL=uL,
    overwrite=True
)

## results
param_names, opt_params, LL = ret

out_params_lkl = open(outfile_lkl, mode="w", encoding="utf-8")

outfile_lkl.write('Log-likelihood: {}\n'.format(-LL))
outfile_lkl.write("Best fit parameters:\n")

for n, p in zip(param_names, opt_params):
    outfile_lkl.write(f"{n}\t{p:.3}\n")
outfile_lkl.close()

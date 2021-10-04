import moments
import demes
import demesdraw
import matplotlib.pyplot as plt
import pickle
import gzip


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
print("Log-likelihood:", -LL)
print("Best fit parameters")
for n, p in zip(param_names, opt_params):
    print(f"{n}\t{p:.3}")





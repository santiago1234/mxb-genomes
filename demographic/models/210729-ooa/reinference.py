import moments
import demes
import demesdraw
import matplotlib.pyplot as plt
import pickle
import gzip


# start with the best model
model_yml = "results/ooa_best_fit_model.yml"
options_params = "ooa-5pops-options.yml"
inferred_model_yml = "results/reinfered-best.yml"

# data
sf = "../../data/210804-Compute-jSFS/data/spectrums/5d-csq-synonymous-spectrum.pkl.gz"
with gzip.open(sf, "rb") as f:
    sf = pickle.load(f)

print(sf.pop_ids)



## project an fold
## We project to a sample size
## such that the MXL population has more samples
n = 60
projection = [n, n, n]
sf = sf.project(projection)
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
    verbose=1,
    maxiter=300,
    uL=uL,
    overwrite=True
)

## results
param_names, opt_params, LL = ret
print("Log-likelihood:", -LL)
print("Best fit parameters")
for n, p in zip(param_names, opt_params):
    print(f"{n}\t{p:.3}")

# Forward in time genetic simulation under inferred demography


## Protocol

I am using the mode [ADMIXTURE-MXL.yml](./ADMIXTURE-MXL.yml),
this model is in generation time see [model-in-generations.py](./model-in-generations.py).
I also round up to integers the MXL start time (16 generations)
and the YRI-AFR pulse in MXL (13 generations).

1. Run the simulations

Running this step with the 30 cores took ~1 month.

```bash
conda activate fwdpy11_18
snakemake -j30  results/simulations/sim-{1..350}-pop.bin
```

NOTE: Follow the order I am indicating.

2. Process the simulations

```
snakemake -j30 results/all-SFSs.csv
```

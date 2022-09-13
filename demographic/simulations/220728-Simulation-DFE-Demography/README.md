# Forward in time genetic simulation under inferred demography


## Protocol

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

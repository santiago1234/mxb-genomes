# LD-Decay-simulated-data
### Santiago G. Medina Muñoz
### 23/05/15

## Overview

The goal of this analysis is to Compute the LD decay curves in the simulated data (forward in time).

## Data

We use the simulated data with [fwdpy11](../../simulations/220728-Simulation-DFE-Demography/README.md)

## Results


```bash
mkdir -p data/recomb_map/
mkdir -p data/vcfs/
mkdir -p data/intervals

snakemake -j1 Snakefile
```


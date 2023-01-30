# PlotSFSExpectedObservedForCategories
### Santiago G. Medina Muñoz
### 23/01/30

## Overview

The goal of this analysis is to Plot the SFS Expected VS Observed across the SNP categories.

The same as [../220204-Expected-vs-Observed-SFS/](../220204-Expected-vs-Observed-SFS/) but for 
all the categories.

## Data

A quick description of the data needed to run this analysis.

- `../../data/220113-ConstructBoostrapedDatasets/data/whole-genome/spectrum-cat_intergenic.pkl.gz` -> 
    The SFS for the whole genome.
- `../220124-InfereModels/results/best-guest-NAT-EXPANSION-intergenic.yml` ->
    The inferred demographic model.

## Results

### Results data

* `data 1`, This data set contains (summary stats etc.)
* `data 2`

### Plots

* `plots/plot1` This plot is ..
* `plots/plot2` This plot is ..

## Conlusions

* Conlusion 1
* Conlusion 2


## Next Steps & Ideas for Future Work

* Next I will do ...


## Reproducibility

To replicate this analysis run:

```bash
snakemake -j1 Snakefile
```


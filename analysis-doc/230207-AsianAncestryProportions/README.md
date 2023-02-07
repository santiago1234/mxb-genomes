# AsianAncestryProportions
### Santiago G. Medina Muñoz
### 23/02/07

## Overview

The goal of this analysis is to Compute the East Asian ancestry
proportions in the Cohorts from ADMIXTURE and Local ancestry.

## Data

See the scripts in [scripts](./scripts) for input datasets.

## Results

### Results data

* `results/ADMIXTURE-ancestry-props.csv`, ancestry proportions from
  ADMIXTURE.
* `results/lai-ancestry-props.csv`, ancestry proportions from local
  ancestry.


## Conlusions

* I detect higher levels from EAS ancestry in MXL than the others with
  ADMIXTURE.
  * MXL ranks top 2 with LAI, PEL is top 1.


## Reproducibility

To replicate this analysis run:


```bash
snakemake -j1 all
```


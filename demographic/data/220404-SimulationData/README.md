## Data for running forward in time genetic simulations

We divided the genome in 1Mb chunks.

NOTE: We wont simulate the whole genome, I will take samples of 1Mb size
only from chromosome 1 and 2.

## Input data

- Chromosome size: [here](../220113-ConstructBoostrapedDatasets/data/human-autosomes.genome)
- Recombination map: See note a the end
- Exon bed intervals: [this script](../211231-mL-coding/data/regions/exones.bed)
- Intronic and intergenic bed intervals: [this script](../211128-compute-mL/00-get-intronic-intergenic-regions.R)



## What do you need to run this pipeline?

- VEP
- bedtools


NOTE: Recombination map.

You can read the recombination map with msprime using

```pytoh
import msprime

rmap = '../../../resources/genetic-maps/chr1.b38.gmap'
rmap = msprime.RateMap.read_hapmap(rmap, position_col=0, map_col=2)
## you can take an slice of the map with 
rmap.slice
```

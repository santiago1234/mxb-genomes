## Data for running forward in time genetic simulations

We divided the genome in 1Mb chunks.

NOTE: We wont simulate the whole genome, I will take samples of 1Mb size
only from.

I will take 200 samples, meaning 200Mb.

## Input data

- Chromosome size: [here](../220113-ConstructBoostrapedDatasets/data/human-autosomes.genome)
- Recombination map: See note a the end
- Intronic and intergenic bed intervals: [this script](../211128-compute-mL/00-get-intronic-intergenic-regions.R)
- Interval for coding regions of the genome: [this script](../220423-Coding-mL-Whole-Genome-test/scripts/get-coding-regions.R)
NOTE: I am ussing the intervals for the coding regions, before i used the intervals for the exonic regions, but
exonic regions may contain non coding exons, for example long non coding rnas.

## What do you need to run this pipeline?

- VEP
- bedtools
- pybedtools

## What I do?

0. Sample regions of the genome to simulate, I sample from the top 3rd regions with higher exon content.
1. Compute mutation rates, here I use the pipeline to compute the mutation rates in coding and 
non coding sequences.
2. Generate bed intervals, in each region, for the exonic and non-coding (intronic and intergenic).
3. Get the recombination map.


NOTE: Recombination map.

You can read the recombination map with msprime using

```pytoh
import msprime

rmap = '../../../resources/genetic-maps/chr1.b38.gmap'
rmap = msprime.RateMap.read_hapmap(rmap, position_col=0, map_col=2)
## you can take an slice of the map with 
rmap.slice
```

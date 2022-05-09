


To get the data I run:

```bash
mkdir test-data
cp ../../../resources/genetic-maps/chr22.b38.gmap test-data
cp ../../data/220404-SimulationData/data/samples/region_*_4.* test-data/
```

## Data

* *test-data/region_region_4.bed* -> The sampled 1Mb genomic region
* *test-data/chr22.b38.gmap* -> Recombination map in chr22
* *test-data/region_exons_4.bed* -> Exonic intervals
* *test-data/region_intronANDinterg_4.bed* -> Non-Exonic intervals (intronic and intergenic)
* *test-data/region_mlcoding_4.csv* -> Scaled mutation rates for missense, synonymous, and loss of function variants
* *test-data/region_mlnoncoding_4.txt* -> Scaled mutation rates for non-coding variants (intronic and intergenic)


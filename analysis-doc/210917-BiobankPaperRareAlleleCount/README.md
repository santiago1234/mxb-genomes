## Analysis of allele counts by NAT ancestry in 1K genomes data 


## Protocol

## Step 1: Liftover bed (MEGA-array regions) to GRCh38

The data I have (1K genomes vcfs) is on GRCh38 the MEGA array regions
are in GRCh37. Then I need to put them in GRCh38.

```
conda activate ggcrossmap
snakemake -j1 data/mega-array-regions/Multi-EthnicGlobal_D1-GRCh38.bed
```



## Mutation Burden Concordance Between the MEGA array and WGS data


NOTE: This analysis is restircted to variants with MAF < 0.05

```bash
snakemake -j30 results/alt-counts.csv
snakemake -j1 clear
```
### Plan

- Count the number of alternative allels in each individual for the categories:
    - synonymous
    * missense
    * LOF

* Do the same in the MEGA regions

* Is there a correlation, which means a strong mutation burden concordance.




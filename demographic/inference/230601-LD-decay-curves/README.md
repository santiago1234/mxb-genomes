# LD decay curves

This analysis is identicall to [this one](../230511-LD-decay-curves/README.md),
but here we use the regions of the geome that are 10Mb in length.

We will focus only on intervals that fall within a single chromosome since the
script takes a vcf from a chromosome as input.

```bash
snakemake -j15 results/ld_stats/{YRI,IBS,CHB,MXB,MXL}-region{1..100}-ld_stats.pkl
```

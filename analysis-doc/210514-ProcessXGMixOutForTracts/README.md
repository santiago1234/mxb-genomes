## Pipeline to process XGMix output for running Tracts


## Step 1: Remove low density windows

Here, I removed the top 22 windows with lowest SNP density. Approx one region per chromosome.

## Step 2: Merge continuous ancestry windows into tracts

## Step 3: To bed like format for running tracts

The data to run tracs is here: *data/3-pops/tracts/PUR/*

To generate the data use:

```
snakemake -j10 all
```

## Step 4: Visualization

Individual karyograms

To generate a particular karyogram run:

```
snakemake -j1 plots/3-pops/karyogram-HG01893.png
```

![karyo](plots/3-pops/karyogram-HG01893.png)

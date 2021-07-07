#!/bin/basg
rm -rf data  results
snakemake -j4 all
snakemake -j1 results/final-sfs.csv
python format-results.py

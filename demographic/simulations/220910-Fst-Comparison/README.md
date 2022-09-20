# Compare Fst values between the data and the model


For each chunk (350) compute the Fst between each pair of
populations.

We will do this for non-coding data.

NOTE: I need to remove mask in the simulation.

For running the pipeline.

```bash
# to generate the fst for data I need
# popgene env
conda activate popgene
snakemake -j10 data_fst


# to generate fst for the simulation I need
# the fwdpy11_18 env
conda activate fwdpy11_18
scnakemake -j10 al

```

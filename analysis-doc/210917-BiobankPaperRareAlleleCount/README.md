## Analysis of allele counts by NAT ancestry in 1K genomes data 


## Protocol

## Step 1: Liftover bed (MEGA-array regions) to GRCh38

The data I have (1K genomes vcfs) is on GRCh38 the MEGA array regions
are in GRCh37. Then I need to put them in GRCh38.

```
conda activate ggcrossmap
snakemake -j1 data/mega-array-regions/Multi-EthnicGlobal_D1-GRCh38.bed
```


## Notes:


Hi Mashal,


This is the updated plot. 


The points correspond to the MXL individual from 1KGP.

Variants were annotated with VEP.

The NAT ancestry proportion was determined with ADMXITURE (k=3), taking the NAT ancestry  component (x-axis).

We polarize for the ancestral allele and removed variants were the ancestral allele was low confidence (lower case, e.g atgc) and now known (“.” Or N).

The y-axis shows the number of derived variants per individual.
-	Homozygous for ancestral allele were counted as 0
-	Heterozygous counted as 1.
-	Homozygous for derived allele as 2.

Then we aggregated (sum) the counts across all variants (22 autosomes) per individual.

Deleterious variants consist of variants with a predicted consequence (in the fitness) of moderate and high: this includes the following variant categories:

	- transcript_ablation
	- splice_acceptor_variant
	- splice_donor_variant
	- stop_gained
	- frameshift_variant
	- stop_lost
	- start_lost
	- transcript_amplificatio
	- inframe_insertion
	- inframe_deletion
	- missense_variant
	- protein_altering_variant



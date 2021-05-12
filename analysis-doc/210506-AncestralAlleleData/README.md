# Ancestral Allele Data


## Overview

Obtain the ancestral allel for the set of variants (1TGP + 50MXB).

Here, I use the [Ancestral Allele](https://github.com/Ensembl/VEP_plugins/blob/release/104/AncestralAllele.pm) vep 
plugin to obtain the data.

The output is a data frame containing the ancestral allel and
the variatn id. Data is saved at: data/aa-chr{chrn}.csv

Example:

| ID              | AA   |
|:----------------|:-----|
| 22:10519389:T:C | T    |
| 22:10519410:G:A | G    |
| 22:10519413:C:G | C    |
| 22:10519438:A:G | A    |
| 22:10519453:T:C | T    |
| 22:10519515:T:C | T    |
| 22:10519530:T:C | T    |
| 22:10519634:T:G | T    |
| 22:10519639:G:C | G    |


## NOTE:

The ancestral allel data generated here, is specific to be used with the input VCF file used to retrieve the data.
If you wish to use other set of variants, you need to re-run the pipeline changin the input vcf in the [Snakefile](Snakefile).
The reason for this, is that the output tables will only contain the ancestral allele the for variants in the input VCF.

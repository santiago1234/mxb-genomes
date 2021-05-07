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

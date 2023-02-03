# ParamsLatexTable
### Santiago G. Medina Muñoz
### 23/02/02

## Overview

The goal of this analysis is to Make a nice table for the inferred parameters to display in the paper.

## Data

A quick description of the data needed to run this analysis.

* `data/descriptions.csv` -> I manually created this file wich contains
  the parameter names and description, to be included in the table.
- `../220124-InfereModels/results/ConfidenceIntervals/OOA-intergenic.tsv`
  File with infered parameters and confidence intervals for the OOA model
- `../220124-InfereModels/results/ConfidenceIntervals/NAT-EXPANSION-intergenic.tsv`
  File with infered parameters and confidence intervals for the NAT model

I generated a copy of these files in the data directory, see the
[Snakefile](./Snakefile).


## Results

### Results data

* A table in latex code with the parameters. 



## Next Steps & Ideas for Future Work

* I will put this table in the paper.
See the notes at the end for some tips of editing
the table with vim.


## Reproducibility

To replicate this analysis run:

```bash
snakemake -j1 Snakefile all
# I run the script in Rstudio and paste
# the output of the last line in the latex supplement
Rscript scripts/make_latex_tabel.R
```


## Notes

- I put the table in the suplement file, I used the following vim commands:
	* Make shorter column names `s/Optimal Value \w*/Value/g` and `s/Standard error \w*/SE/g`
- Fix columns so they are interpreted in math mode:
	* `s/\\\$/$/g`
	* `s/\\_/_`
	* `s/\\{/{/g`
	* `s/\\}/}/g`
	* `s/textbackslash{}//g`
	* `s/e-\(\d*\)/\\times 10\^\{-\1\}/g` this fixes a better exponential notation.
- After this i made some manual changes.


# Demographic Modeling of Admixed Latin American Populations from Whole Genomes

This repository contains the code and data used in our study, "Demographic
Modeling of Admixed Latin American Populations from Whole Genomes". In this
study, we developed novel demographic models that capture the complex
evolutionary history of Latin American populations, which have been shaped by
both recent admixture and deeper-in-time demographic events. Our models improve
upon existing approaches to infer patterns of genetic variation in admixed
populations, and can be a valuable resource for improving genomic studies and
prediction in these populations.

## Data

We used high-coverage whole genome data from Indigenous American ancestries in
present-day Mexico, as well as existing genomes from across Latin America, to
infer multiple demographic models. The data used in our study is available from
public repositories and is described in the paper.


## Code

The code used to analyze the data and infer demographic models is provided in
this repository. We used a combination of analyses of allele frequencies and
ancestry tract length distributions to develop our models. The code is written
in Python and makes use of several open-source libraries, including `demes`,
`msprime`, `moments`, `tracrs`, and `scikit-allel`.

## Demographic Models

We provide the inferred demographic models in the models directory, in the demes format. The models include:

- [m1-out-of-africa.yml](models/m1-out-of-africa.yml): Four populations out of Africa.
- [m2-Mexico-admixture.yml](models/m2-Mexico-admixture.yml): Admixture in Mexico.
- [m3-Colombia-admixture.yml](models/m3-Colombia-admixture.yml): Admixture in Colombia.
- [m4-Peru-admixture.yml](models/m4-Peru-admixture.yml): Admixture in Peru.
- [m5-PuertoRico-admixture.yml](models/m5-PuertoRico-admixture.yml): Admixture in Puerto Rico.
- [m6-All-admixture.yml](models/m6-All-admixture.yml): Admixture in Latin American populations.

## Usage

To use this Models, you can install the demes library and load the models using the demes.load function.

```Python
import demes

model = demes.load("models/m1-out-of-africa.yml")
```

## Citation

If you use any part of this code or the inferred models in your work, please cite our paper:

```
[Add citation here]
```

Please let us know if you have any questions or issues with the code or data provided in this repository.



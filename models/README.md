# Inferred demographic models.


All the models presented here are in the
[demes](https://academic.oup.com/genetics/article/222/3/iyac131/6730747)
format, a standard format for demographic models.


The inferences from allele frequencies were obtained from the intronic site
frequency spectrum. We obtained similar inferences with intergenic and
synonymous sites.


## Inferred Models

* Model 1 `m1-out-of-africa.yml`: Four populations out of Africa.
* Model 2 `m2-Mexico-admixture.yml`: Admixture in Mexico.
* Model 3 `m3-Colombia-admixture.yml`: Admixture in Colombia.
* Model 4 `m4-Peru-admixture.yml`: Admixture in Peru.
* Model 5 `m5-PuertoRico-admixture.yml`: Admixture in Puerto Rico.
* Model 6 `m6-All-admixture.yml`: Admixture all populations.

See [Models-Visualization](Models-Visualization.ipynb) for a visualization of these models.

***

NOTES (for myself):

To get the model:

```bash
cp ../demographic/inference/220225-Combine-Inferences/ADMIXTURE.yml m6-All-admixture.yml
```

Then I added some manual edits according to each model.

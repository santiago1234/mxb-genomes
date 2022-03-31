# Inference of admixture history


Here, I infere the admixtur history for the populations: MXL, NAT, CLM, and PEL.

I tested the following set of models:

- ppx_xxp
- ppx_xxp_pxx
- ccx_xxp
- ppx_ccx_xxp
- ppp_pxp  (This model is directed for PUR). NOTE: This model wont
be shown. It does not captures well the historical records.
- ppc. This is the model that i decided to use for PUR.

I run some tests before in [analysis-doc](`../../../analysis-doc), but here is the final inference.

The ancestry order is always EUR, NAT, and MXL. But in PUR it is AFR, EUR, and NAT.

I have put these models in a module: [tractsmodels](tractsmodels/).


**NOTES:**

+ When I run the inference with snakemake I get an error. I think the reason is
that tracts messages are interpreted as errors.
+ To solve that I am saving the python calls to a file `runall.sh` and then passing
that to the unix `parallel` commmand.


```bash
snakemake -j52 infere_all
# this writes the python calls to the runall.sh file
cat runall.sh |parallel -j22
```

**IMPORTANT NOTE:** Check the order of the labels for the models, see [this script](tractsmodels/utils.py),
where the labels order is defined.

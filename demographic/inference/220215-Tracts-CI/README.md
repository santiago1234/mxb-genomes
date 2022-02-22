# Inference of admixture history


Here, I infere the admixtur history for the populations: MXL, NAT, CLM, and PEL.

I tested the following set of models:

- ppx_xxp
- ppx_xxp_pxx
* ccx_xxp
* ppx_ccx_xxp

I run some tests before in [analysis-doc](`../../../analysis-doc), but here is the final inference.

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

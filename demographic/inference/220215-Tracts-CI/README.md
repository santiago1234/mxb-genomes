# Estimate Tracts confidence intervals




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

# Construct bootstrapped datasets of the SFS and associated scaled mutation rates

I dive the genome in 384 chunks each spaning ~6Mb (this taking into account the mask).


## jSFS


I compute a 4 dimensional jSFS. With the populations IBS, CHB, MXB, and YRI. The output spectrum
is **not polarized**.

See [this script](scripts/define-pops.py) to see how I defined the samples.


* I removed the CpGs mutations.

## Mutation rates


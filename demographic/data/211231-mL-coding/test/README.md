# Test

I run the test in chr22.

# Test 1: Non-overlapping variants.

First, I check that the are not repeated variants, in the universe of all possible biallelic variants. In case,
o repeated variants, we could count them twice for mL calculation, which we do not want.

```bash
python test-not-overlaping.py
[output]:
    3    2658201
```


What that means? 
The number **3** says that each position is represented three times. The number 3 is because we have 3 alternative alleles for each position.
The numer 3 is the only number that I expect to observe. The test is sucessful.


# Test 2: Variants observed in the genotypes are a subset of all possible variants (universe).

* I took all variants, in the genotypes, that are missense in chr22.
* Then I generete a list with all the possible variants (universe) that were annotated as missense.

The test consists in checking that the genotypes (missense) variants are a subset of the universe (missense) variants.

```python
python test-vcf-variants-is-a-subset.py
[output]:
    Test failed at position: 33305584
    Test failed at position: 38298367
    Test failed at position: 38298433
    Test failed at position: 40491992
    Test failed at position: 40491993
    Test failed at position: 40492086
    Test failed at position: 40492104
    Test failed at position: 40492152
    Test failed at position: 41507468
    Test failed at position: 46390911
    Test failed at position: 46390998
    Test failed at position: 46391009
    Test failed at position: 46391054
    Test failed at position: 46391056
    Test failed at position: 50698183
    Test failed at position: 50698218
    Test failed at position: 50705397
```

Only 16 variants were not found in the universe list. There are 4815 missense variants in chr22. 
So 99.667% of the genotype variants were contained in the universe list. Since this is most of the 
variants I guess this is a passed test.

I looked more into this, I generated the following bedfile, representing the not found variants:

```
22	33305583	33305585
22	38298366	38298368
22	38298432	38298434
22	40491991	40491993
22	40491992	40491994
22	40492085	40492087
22	40492103	40492105
22	40492151	40492153
22	41507467	41507469
22	46390910	46390912
22	46390997	46390999
22	46391008	46391010
22	46391053	46391055
22	46391055	46391057
22	50698182	50698184
22	50698217	50698219
22	50705396	50705398
```

I took the intersection of this file with the exons file (this after masking).
The intersection was empty which indicates that these variants were not in my starting
exons list. This is a good indication.



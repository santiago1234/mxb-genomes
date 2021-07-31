# Test: Ancestral genome matches VCF files

When I [ran relate](../../analysis-doc/210723-GeneGenealogies), during the polzarization step,
relate drops ~1/2 of the SNPs.
Also, the estimated population sizes look not as what we would expect. Here, I will run some analyses
to make sure the ancestral genome is aligned to the VCF data.

I obtained the ancestral genome from [here](http://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/).

The script [test-ancestral-genome-match.py](test-ancestral-genome-match.py) generates the following stats table:


| case          |      n |   percentage |
|:--------------|-------:|-------------:|
| ALT_is_ANC    |  25140 |      7.14477 |
| REF_is_ANC    | 304360 |     86.4988  |
| ANC_not_KNWON |  22366 |      6.3564  |


This result looks fine, so the ancestral genome is aligned to the VCF file(s).

## Conclusion

I think realtes expects one and only one ancestral sequence in the fasta file. (Currently) I am
passing a fasta file with the ancestral sequences for all the genome (all chromosomes).

I was generating this ancestral genome with:


```bash
wget ftp://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz
tar xfz homo_sapiens_ancestor_GRCh38.tar.gz
cat homo_sapiens_ancestor_GRCh38/*.fa | bgzip -c > homo_sapiens_ancestor_GRCh38.fa.gz
rm -rf homo_sapiens_ancestor_GRCh38/ homo_sapiens_ancestor_GRCh38.tar.gz
mv homo_sapiens_ancestor_GRCh38.fa.gz 210719-ancestral-genome/
```

I conclusion, I should past to Relate a fasta file with only the ancestral sequence for one chromosome (one record).


There are a few duplicated variants.
I removed this variants with

```bash
bcftools norm --remove-duplicates variant_effect_output.vcf.gz  -Oz -o clean.vcf.gz
```

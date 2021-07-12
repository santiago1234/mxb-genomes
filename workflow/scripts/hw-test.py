"""
Compute the hw for the variants in a vcf file.
The output is a csv file with columns:
    + ID -> the variant ID
    + chi2 -> chi square value
    + p -> the p value

usage:
    python hw.py <vcf_file> <outfile.csv>
NOTE:
    varinats were hw cannot be computed will
    not be in the output table. Sometimes HW
    cannot be computed because there is no
    variantion in the SNP (e.g. all 0/0 or 1/1)
"""

import pandas as pd
import allel
import sys
sys.path.append("./")
from mxbgenomes.stats import hw_test

vcf_file = sys.argv[1]
outfile = sys.argv[2]

vcf = allel.read_vcf(vcf_file)
ga = allel.GenotypeArray(vcf['calldata/GT'])
chi, p = hw_test(ga)


vars_id = vcf['variants/ID']
d = {
    'ID': vars_id,
    'chi2': chi,
    'p': p
}

col_order = ['ID', 'chi2', 'p']
(
    pd.DataFrame(d)
    .dropna()
    [col_order]
    .to_csv(outfile, index=False)
)


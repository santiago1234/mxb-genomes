from cyvcf2 import VCF
from Bio import SeqIO
import pandas as pd

path_to_vcf = "../../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr22-snps-vep-mask-HW-GRCh38.vcf.gz"
path_to_ancestral_genome = "../../tmp/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_22.fa"

# Read the ancestral sequence
record = SeqIO.read(path_to_ancestral_genome, "fasta")

# I create this list that will hold the data to later
# create a pandas frame
d = []
for variant in VCF(path_to_vcf):
    vd = [variant.ID, variant.REF, variant.ALT[0]]  # I am assuming the input are only biallelic variants
    # get the ancestral allel
    ancestral_allel = record.seq[variant.start]
    vd.append(ancestral_allel)
    d.append(vd)


result = pd.DataFrame(d, columns=['ID', 'REF', 'ALT', 'ANC'])

# make the ancestral upper-case, from the readme file the lower
# case means that there is less confidence
result['ANC'] = result.ANC.str.upper()


# check cases
alt_is_anc = (result.ALT == result.ANC).sum()
ref_is_anc = (result.REF == result.ANC).sum()
# By not know I mean that the ancestral allel is not
# equal to the REF or the ALT.
anc_not_known = result.shape[0] - alt_is_anc - ref_is_anc

# make output table with summary

d = {
    'ALT_is_ANC': [alt_is_anc],
    'REF_is_ANC': [ref_is_anc],
    'ANC_not_KNWON': [anc_not_known]
}
d = pd.DataFrame(d).melt(var_name='case', value_name='n')
d['percentage'] = 100 * d['n'] / d['n'].sum()

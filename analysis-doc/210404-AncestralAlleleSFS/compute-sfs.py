"""
Compute unfolded SFS from VCF file:
usage:
	python <VCF> <AA> <output> <statsfile>
Note:
	This script migth only work inside this dir
	due to relative paths used.
"""
import sys
sys.path.append("../..")
from mxbgenomes import sfs
from mxbgenomes.utils import load_populations_info

vcf_file = sys.argv[1]
aa_file = sys.argv[2]
outfile = sys.argv[4]
statsfile = sys.argv[4]

popinfo = load_populations_info("../../")

subpops = (
    popinfo
    .groupby('Subpopulation')
    .apply(lambda x: x.Samplename.to_list())
    .to_dict()
)

spectrum, stats = sfs.sfs_unfolded(vcf_file, aa_file, subpops=subpops, project_haplod_size=50)
spectrum, stats = sfs.sfs_to_frame(spectrum)
spectrum.to_csv(outfile, index=False)
stats.to_csv(statsfile, index=False)

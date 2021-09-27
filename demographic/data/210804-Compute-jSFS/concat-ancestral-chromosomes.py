from Bio import SeqIO
import os
import glob



def load_chr(n):
    """
    load the ancestal chromosome into a seqrecord
    n: int 2-22
    """
    path_to_ancestral = "../../../../gene-genealogies-mxb/resources/210719-ancestral-genome/homo_sapiens_ancestor_GRCh38/"
    fname = "homo_sapiens_ancestor_" + str(n) + ".fa"
    fname = os.path.join(path_to_ancestral, fname)
    print('loading: {}'.format(fname))
    record = SeqIO.read(fname, format="fasta")
    record.id = str(n)
    return record


records = [load_chr(i) for i in range(1, 23)]
outdir = 'data/ancestral-genome'
if not os.path.exists(outdir):
    os.makedirs(outdir)

outfile = os.path.join(outdir, "ancestral-genome-autosomes.fasta")
SeqIO.write(records, outfile, "fasta")






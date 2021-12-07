"""
script to count 3 mers across a set of sequences.
The results is the aggregated total kmer count.

usage:
    python count-context.py <ses.fasta> <outputTable.csv> <cores>
"""
import sys
import multiprocessing
import functools
from collections import Counter
import pandas as pd
from Bio import SeqIO
from countmers import count_3mers_in_record


def get_kmers(aseq, k):
    return [aseq[i: i+k] for i in range(len(aseq) - k + 1)]


def count_kmers(aseq, k):
    return Counter(get_kmers(aseq, k))


def count_3mers_in_records(records_list):
    """
    Count the 3mers in a list of records
    """
    print('processing records >>>>')
    count3mers = lambda x: count_kmers(str(x.seq).upper(), 3)
    cuentas = [count3mers(x) for x in records_list]
    return functools.reduce(lambda a, b: a + b, cuentas)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def run(seqs_fasta, output_tab, cores):
    """
    seqs_fasta: str, path to file with sequences to count 3-mers
    """
    records = SeqIO.parse(seqs_fasta, "fasta")
    records = list(records)
    # make sublists for parallel processing
    records = list(chunks(records, 1000))
    print('>>> COMPUTING <<<')
    pool = multiprocessing.Pool(processes=cores)
    kcounts = pool.map(count_3mers_in_records, records)
    kcounts = functools.reduce(lambda a, b: a + b, kcounts)
    # make a data fram
    kcounts = (
        pd.DataFrame
        .from_dict(kcounts, orient='index')
        .reset_index()
        .rename(columns={'index': 'kmer', 0: 'count'})
    )
    kcounts.to_csv(output_tab, index=False)
    pass


if __name__ ==  '__main__':
    seqs_fasta, output_tab, cores = sys.argv[1:]
    cores = int(cores)
    run(seqs_fasta, output_tab, cores)

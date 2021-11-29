import multiprocessing
from collections import Counter
from Bio import SeqIO
from countmers import count_3mers_in_record


aseq = "ATGCCGTA"
cores = 4

if __name__ ==  '__main__':
    records = SeqIO.parse("data/regions/introns.fasta", "fasta")
    records = list(records)
    print('>>> COMPUTING <<<')
    pool = multiprocessing.Pool(processes=cores)
    kcounts = pool.map(count_3mers_in_record, records)

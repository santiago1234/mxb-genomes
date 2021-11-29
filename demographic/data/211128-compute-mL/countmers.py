from collections import Counter

def get_kmers(aseq, k):
    return [aseq[i: i+k] for i in range(len(aseq) - k + 1)] 

def count_kmers(aseq, k):
    return Counter(get_kmers(aseq, k))

def count_3mers_in_record(record):
    aseq = str(record.seq)
    print(f'Process {record.id} >>>')
    return count_kmers(aseq.upper(), 3)

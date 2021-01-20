#! /usr/bin/env python3

import sys
from Bio import SeqIO


file = sys.argv[1]
with open(file, 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(f'{record.id} {len(record.seq)}')

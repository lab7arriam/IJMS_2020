#! /usr/bin/env python3

import os
import sys
from Bio import SeqIO


dir = sys.argv[1]

for file in os.listdir(dir):
    with open(os.path.join(dir, file), 'r') as handle:
        counter = 1
        for record in SeqIO.parse(handle, 'fasta'):
            print(f'{file}({counter}): {len(record.seq)}')
            counter += 1

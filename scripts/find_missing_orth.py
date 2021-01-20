#! /usr/bin/env python3

import os
from Bio import SeqIO

proteomes  = os.path.join(os.path.expanduser('~'), 'bacillus_phylogeny_2020', 'flagellin_phylogeny', 'proteomes')
prefix = 'GCA_000161535'
protein = 'bthur0004_49390'
for file in os.listdir(proteomes):
    if file.find(prefix) > -1:
        with open(os.path.join(proteomes, file), 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if record.id.find(protein) > -1:
                    with open(f'{prefix}.missing_nprb.faa', 'w') as handle2:
                        SeqIO.write(record, handle2, 'fasta')

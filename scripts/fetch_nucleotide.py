#! /usr/bin/env python3

import click
import re
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

@click.command()
@click.option('-n', default='.')
@click.option('-d', default=None)
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path(exists=False))


def find_nucleotides(input, output, n, d):

    '''Create an output directory'''

    if not os.path.exists(output):
        os.makedirs(output)

    '''Parse Roary cluster data'''

    roary = pd.read_csv(d, header=0, sep=',')
    rcols  = [x for x in roary.columns.values if x.find('GCA') > -1]
    end = len(roary.index)
    clusdict = dict()

    for rcol in rcols:
        for i in range(0, end):
            orth_string = roary[rcol][i]
            if not pd.isnull(orth_string):
                for orth in orth_string.split('\t'):
                    clusdict[orth] = roary.iloc[i,0]


    '''Parse an HMMer/BLAST output and extract all sequence IDs,
    then gather all nucleotide sequence files and write them to the output directory'''

    files = os.listdir(input)
    nuc_seqs = os.listdir(n)

    for file in files:
        prefix = file.split('.')[0]
        file_list = set()
        output_list = []
        with open(os.path.join(input, file), 'r') as handle:
            for line in handle.readlines():
                if line.find('#') == -1:
                    gene = line.split(' ')[0]
                    file_list.add(gene)
        for id in file_list:
            key = clusdict[id]
            for nuc_id in nuc_seqs:
                if nuc_id.find(key) > -1: 
                    with open(os.path.join(n, nuc_id), 'r') as handle:
                        records = SeqIO.parse(handle, 'fasta')
                        for record in records:
                            if record.id == id:
                                record.seq = Seq(str(record.seq).replace('-', ''))
                                record.description = f'cluster={key}'
                                output_list.append(record)

        output_name = f'{prefix}.fasta'
        with open(os.path.join(output, output_name), 'w') as handle:
            SeqIO.write(output_list, handle, 'fasta')


if __name__ == '__main__':
    find_nucleotides()

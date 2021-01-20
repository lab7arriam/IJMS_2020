#! /usr/bin/env python3

import click
from Bio import SeqIO

@click.command()
@click.option('-i', type=str)
@click.option('-n', type=str)
@click.option('-o', type=str)


def extract_proteins(i, n, o):

    with open(n, 'r') as handle:
        record_list = list(map(lambda x: x.strip('\n'), handle.readlines()))

    print(record_list)
    output_list = []

    with open(i, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in record_list:
                output_list.append(record)

    with open(o, 'w') as handle:
        SeqIO.write(output_list, handle, 'fasta')

if __name__ == '__main__':
    extract_proteins()

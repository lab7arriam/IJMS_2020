#! /usr/bin/env python3

import click
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
Get the longest/shortest sequence from the list of orthologs
"""
def extractByLength(input, query, output):

    with open(input, 'r') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
        if len(records) > 0:
            if query == 'long':
                record = max(records, key=lambda x: len(x.seq))
            elif query == 'short':
                record = min(records, key=lambda x: len(x.seq))
        else:
            record = SeqRecord(id = 'None', seq=Seq(''))
        with open(output, 'w') as handle:
            SeqIO.write(record, handle, 'fasta')


@click.command()
@click.option('--input', '-i', help='An input file containing orthologous sequences in FASTA format', type=click.File('rb'))
@click.option('--query', '-q', help='Define whether the longest/shrortest sequence should be extracted',  type=click.Choice(['long', 'short'], case_sensitive=False))
@click.option('--output', '-o', help='A path to the output file', type=click.File('wb'))

def extract(input:str, query:str, output:str):

    if not os.path.exists(output):
        ps.makedirs(output)

    extractByLength(input, query, output)


if __name__ == '__fasta__':
    extract()

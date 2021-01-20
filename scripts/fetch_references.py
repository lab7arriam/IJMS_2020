#! /usr/bin/env python3

import click
import os
import sys
import pandas as pd
from Bio import SeqIO

@click.command()
@click.option('--email', '-e', help='An email address to introduce to NCBI with', type=str)
@click.option('--input', '-i', help='A path to the input file', type=str)
@click.option('--output', '-o', help='A path to the output directory')

def fetchReferences(email, input, output):

    """
    Reads a headless accession file (one accession in a row), creates a temporary table file and sends to fetchSequences script
    """

    if not os.path.exists(output):
        os.makedirs(output)

    with open(input, 'r') as handle:
        references = pd.read_csv(handle, header=None, sep='\t')

    for i in range(0, references[0].count()):
        tmp_df = pd.DataFrame({'Name':[references.iloc[i,0]]})
        with open('tmp_df.tsv', 'w') as handle:
            tmp_df.to_csv(handle, header=True, index=False, sep='\t')

        filename = f'spot{references.iloc[i,1]}_{references.iloc[i,0]}.fasta'
        output_path = os.path.join(output, filename)

        sys.path.append('/home/yura/bacillus_phylogeny_2020/DIGEs/')
        from fetch_sequences import fetchSequences
        fetchSequences(email=email, input='tmp_df.tsv', output=output_path)

        os.remove('tmp_df.tsv')


if __name__ == '__main__':
    fetchReferences()

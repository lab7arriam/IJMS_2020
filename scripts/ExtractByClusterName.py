#! /usr/bin/env python3

import click
import os
import pandas as pd
from Bio import SeqIO

"""
Extract nucleotide sequences by accessions of  Roary clusters representative sequences.
"""

def extractAccessions(input:str, accol:int, namecol:int)->dict:

    with open(input, 'r') as handle:
        clust_dict = dict()
        df = pd.read_csv(handle, header=0, sep='\t')
        for i in range(0, len(df)):
            clust_dict[df.iloc[i, accol].split('.')[0]] = df.iloc[i, namecol]

        return clust_dict


def extractSequence(input:'generator', acc:str, cluster_dict:dict)->'Bio.SeqRecord.SeqRecord':

    for record in input:
        if record.description.split('=')[-1] == cluster_dict[acc]:
            return record


def globalExtract(indir:str, outdir:str, clust_dict:dict):

    files = os.listdir(indir)
    for file in files:
        prefix = file.split('.')[0]
        with open(os.path.join(indir, file), 'r') as handle:
            sequences = SeqIO.parse(handle, 'fasta')
            record = extractSequence(sequences, prefix, clust_dict)

            filename = f'{prefix}.fasta'
            with open(os.path.join(outdir, file), 'w') as handle:
                SeqIO.write(record, handle, 'fasta')


@click.command()
@click.option('--input', '-i', help='A path to the input directory', type=click.Path(exists=True))
@click.option('--output', '-o', help='A path to the output directory', type=click.Path(exists=False))
@click.option('--clusters', '-c', help='A tabular-separated table containing cluster-representing sequences listed in an accession-specific manner', type=str)
@click.option('--accessions', '-a', help='A number defining the column containing genome accessions (default: 0)', default=0, type=int)
@click.option('--names', '-n', help="A number defining the column containing representative sequences' accessions", type=int)


def launch(input, output, clusters, accessions, names):

    if not os.path.exists(output):
        os.makedirs(output)

    clust_dict = extractAccessions(clusters, accessions, names)
    globalExtract(input, output, clust_dict)


if __name__ == '__main__':
    launch()

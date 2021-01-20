#! /usr/bin/env python3

import click
import os
import regex
import pandas as pd
from Bio import SeqIO

"""
A crutchy script for adding the longest amplicon cluster to the existing amplicon cluster pivot
"""

def extractCluster(input:str)->str:

    with open(input, 'r') as handle:
        record = SeqIO.read(handle, 'fasta')
        clusregex = regex.compile('cluster=[\.A-Z0-9_]+', flags=regex.IGNORECASE)
        cluster = regex.findall(clusregex, record.description)

        return cluster[0]


def globalExtractCluster(input:str):

    files = os.listdir(input)
    cluster_dict = dict()
    for file in files:
        nameregex = regex.compile('GCA_[0-9]+')
        prefix = regex.findall(nameregex, file)[0]
        cluster = extractCluster(os.path.join(input, file))
        cluster_dict[prefix] = cluster.replace('cluster=', '')

    cluster_df = pd.DataFrame({'Accession': list(cluster_dict.keys()), 'Cluster': list(cluster_dict.values())})

    return cluster_df

def addColumn(main:pd.core.frame.DataFrame, clusters:pd.core.frame.DataFrame):

    sorter = main[main.columns.values[0]]

#    print(clusters)

    clusters.Accession = clusters.Accession.astype('category')
    clusters.Accession.cat.set_categories(sorter, inplace=True)
    clusters.sort_values(['Accession'], inplace=True)

#    print(clusters)

    main.insert(main.shape[1], 'PCR_Ref_Cluster', clusters['Cluster'].to_list(), True)

    return main


@click.command()
@click.option('--table', '-t', type=click.File())
@click.option('--input', '-i', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(exists=False))


def execute(table, input, output):

    clusters = globalExtractCluster(input)
    table = pd.read_csv(table, sep='\t', header=0)
    output_df = addColumn(table, clusters)

    output_df.to_csv(output, sep='\t', header=True, index=False)

if __name__ == '__main__':
    execute()

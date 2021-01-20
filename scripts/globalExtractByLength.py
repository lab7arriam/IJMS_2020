#! /usr/bin/env python3

import click
import os
from ExtractByLength import extractByLength

"""
extract longest/shortest orthologs from a load of FASTA sequence files stored in on directory
"""

def extractGlobal(input, query, output):

    if not os.path.exists(output):
        os.makedirs(output)

    for file in os.listdir(input):
        outname = f'{file.replace(".fasta", "").replace(".fa", "")}_{query}est.fasta'
        print(file)
        extractByLength(input=os.path.join(input, file), query=query, output=os.path.join(output, outname))


@click.command()
@click.option('--input', '-i', help='A path to the directory containing input files', type=click.Path(exists=True))
@click.option('--query', '-q', help='Define whether the longest/shortest isoform should be extracted from each file', type=click.Choice(['long', 'short'], case_sensitive=False))
@click.option('--output', '-o', help='A path to the output directory (creates one if non-existent)', type=click.Path(exists=False))

def extract(input:str, query:str, output:str):

    extractGlobal(input, query, output)


if __name__ == '__main__':
    extract()

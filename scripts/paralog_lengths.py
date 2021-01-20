#! /usr/bin/env python3

import click
import os
import pandas as pd
from Bio import SeqIO

@click.command()
@click.argument('input_dir')
@click.argument('output_file')


def check_lengths(input_dir, output_file):

    names, lengths = [], []

    for file in os.listdir(input_dir):
        alias = '_'.join(file.split('_')[0:2])
        names.append(alias)
        length_list = []
        with open(os.path.join(input_dir, file), 'r') as handle:
            for entry in SeqIO.parse(handle, 'fasta'):
                length_list.append(str(len(str(entry.seq))))
        lengths.append(';'.join(length_list))

    output = pd.DataFrame({'Accession': names, 'Default_length': lengths})
    output.to_csv(output_file, sep='\t', header=True, index=False)

if __name__ == '__main__':
    check_lengths()

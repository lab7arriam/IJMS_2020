#! /usr/bin/env python3

import click
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

@click.command()
@click.option('--blast_dir', '-bd', type=str)
@click.option('--prot_dir', '-pd', type=str)
@click.option('--gene-name', '-g', type=str, default='hag')
@click.option('--output', '-o', type=str, default='.')

def extractProteins(blast_dir, prot_dir, gene_name, output):

    if not os.path.exists(output):
        os.makedirs(output)

    blast_results = os.listdir(blast_dir)
    for file in blast_results:
        gene_list = []
        output_list = []
        prefix = file.split('.')[0]
        with open(os.path.join(*[blast_dir, file]), 'r') as handle:
            for line in handle.readlines():
                gene_list.append(line.split('\t')[0])
#            prot_file = f'{filename.rstrip(".proteome")}_{gene_name}.faa'
            for prot_file in os.listdir(prot_dir):
                if prot_file.find(prefix) > -1:
                    with open(os.path.join(prot_dir, prot_file), 'r') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            for entry in gene_list:
                                if record.id == entry:
                                    output_list.append(record)
            output_file = f'{prefix}.{gene_name}.fasta'
            with open(os.path.join(output, output_file), 'w') as handle:
                SeqIO.write(output_list, handle, 'fasta')


if __name__ == '__main__':
    extractProteins()

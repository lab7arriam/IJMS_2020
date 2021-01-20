#! /usr/bin/env python3

import os
import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@click.command()
@click.option('--input', '-i', help='Path to input directory', type=click.Path(exists=True))
@click.option('--output', '-o', help='Path to output directory (if non-existent, creates an output directory)', type=click.Path(exists=False))

def generate_fasta_files(input:str, output:str):

    if not os.path.exists(output):
        os.makedirs(output)

    files = os.listdir(input)
    for file in files:
        with open(os.path.join(input, file), 'r') as handle:
            counter = 0
            record_list = []
            prefix = f'spot{file.split(".")[0]}'
            for line in handle.readlines():
                line = line.strip('\n')
                if line.find('(') > -1:
                    line1 = line.split(')')[1].split('(')[0]
                    line2 = line.replace('(', '').replace(')', '')
                    for num, seq in enumerate([line1, line2]):
                        record = SeqRecord(seq=Seq(seq), id=f'{prefix}_seq{counter}_var{num+1}')
                        record_list.append(record)
                else:
                    record = SeqRecord(seq=Seq(seq), id=f'{prefix}_seq{counter}')
                    record_list.append(record)
                counter += 1
            filename = f'{prefix}.fasta'
            with open(os.path.join(output, filename), 'w') as handle:
                SeqIO.write(record_list, handle, 'fasta')


if __name__ == '__main__':
    generate_fasta_files()

#! /usr/bin/env python3

import click
import regex as re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_primers(primers):

    '''Step 3: parse a primer file '''

    output = []
    with open(primers, 'r') as handle:
        for line in handle.readlines():
            line = line.strip('\n')
            local_output = line.split('\t')
            output.append(local_output)

    return output


def generate_templates(primer1, primer2):

    '''Step 1: create a nucleotide alphabet'''

    alphabet = {'W': '[AT]', 'S': '[CG]', 'K': '[GT]', 'R': '[AG]',
                    'Y': '[CT]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]',
                    'V': '[ACG]', 'N': '[ACGT]'}

    forward = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T'}
    reverse = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    forward_al = forward.copy()
    forward_al.update(alphabet)

    reverse_al = reverse.copy()
    reverse_al.update(alphabet)

    '''Step 2: generate primers'''

    primer1 = ''.join(list(map(lambda x: forward_al[x], primer1)))
    primer2 = ''.join(list(map(lambda x: reverse_al[x], primer2[::-1])))

    return [primer1, primer2]


@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path(exists=False))
@click.option('-p', help='Tab-delimited primer file', type=str, default=None)
@click.option('-e', help='Allowed number for primer pair mismatches', type=str, default=4)


def PCR_emulation(input, output, e, p):

    '''Get primers'''
    primers = get_primers(p)

    '''Step 4: check if the output directory exists and create one if missing'''
    if not os.path.exists(output):
        os.makedirs(output)

    '''Step 4: parse sequence files and perform PCR emulation'''
    files = os.listdir(input)
    for file in files:
        prefix = file.split('.')[0]
        output_name = f'{prefix}_PCR.fasta'
        output_list = []
        with open(os.path.join(input, file), 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):

                '''Step 5: generate PCR primers'''

                cprimer = 1
                for primer_pair in primers:
                    primer1, primer2 = generate_templates(*primer_pair)
#                    print(f'Forward: {primer_pair[0]}/{primer1}')
#                    print(f'Reverse: {primer_pair[1]}/{primer2}')
#                    pattern = fr'({primer1}.*?{primer2})' + r'{' + fr'e<={e}' + r'}'
#                    pattern = fr'{primer2}'
                    pattern = fr'{primer1}[ACGT]*{primer2}'
                    regex = re.compile(pattern, re.IGNORECASE)
                    PCR_results = regex.findall(str(record.seq), overlapped=True)

                    '''Step 6: emulate PCR'''

                    cproduct = 1
                    for sequence in PCR_results:
#                        print(len(record.seq), len(sequence))
                        id = f'{record.id}.primer{cprimer}.product{cproduct}'
#                        print(record.description)
                        description = f'{record.description.split(" ")[1]} primer_pair=BthahF{cprimer}/R{cprimer} forward={primer_pair[0]} reverse={primer_pair[1]}'
                        PCR_product = SeqRecord(seq=Seq(sequence), id=id, description=description)
                        output_list.append(PCR_product)
                        cproduct += 1
                    cprimer += 1

        with open(os.path.join(output, output_name), 'w') as handle:
            SeqIO.write(output_list, handle, 'fasta')


if __name__ == '__main__':

#    primers = [['AGTACATGCGCCAAAACCAAG', 'GTTTGCTTGAGAAAGCATGCT'],
#               ['GGGGTTCTTAATCATGAGAA', 'TAACTCAAATGGCTTATTGT'],
#               ['AAYATTAAYAGCATGCGTAC', 'TTTGTGGWGTTTGGTTWGCT']]

    PCR_emulation()

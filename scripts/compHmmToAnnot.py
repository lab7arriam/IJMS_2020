#! /usr/bin/env python3

import click
import math
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

'''

Parse a pivot table containing Roary-deduced orthologs, estimate novel orthologs found by HMMER and extract them to the specified folder

'''


@click.command()
@click.option('-h', type=click.Path(exists=True), help='A directory containing HMMER outputs')
@click.option('-r', type=click.Path(exists=True), help='A Roary outputs pivot table')
@click.option('-o', type=click.Path(exists=False), help='An output folder')

def compare(h, r, o):

    if not os.path.exists(o):
        os.makedirs(o)

#    hmm = pd.read_csv(h, header=0, sep='\t')
    roary = pd.read_csv(r, header=0, sep=',')

    rcols = list(roary.columns.values)
    rind = [i for i in range(len(rcols)) if rcols[i].find('GCA') > -1]

    missing_dict = defaultdict(list)

    for i in range(rind[0], rind[-1]+1):
        roary_prots = roary.iloc[:,i].to_list()
        acc = rcols[i].split('.')[0]
#        print(acc, roary_prots)
        missing_list = []
        for name in os.listdir(h):
            if name.find(acc) > -1:
                with open(os.path.join(h, name), 'r') as handle:
                    for entry in SeqIO.parse(handle, 'fasta'):
                        if entry.id not in roary_prots:
                            missing_dict[acc].append(entry.id)
                            missing_list.append(entry)
                    if len(missing_list) > 0:
                        filename = f'{name.split(".")[0]}.missing'
                        with open(os.path.join(o, filename), 'w') as handle:
                            SeqIO.write(missing_list, handle, 'fasta')
        if len(missing_dict[acc]) == 0:
            missing_dict[acc].append('NA')


    nomiss = list(map(lambda x: len(x), missing_dict.values()))
#    accmiss = list(map(lambda x: ';'.join(x), missing_dict.values()))
    accmiss = [';'.join(x) if x != 'NA' else x for x in missing_dict.values()]
#    print(len(rcols[rind[0]:rind[-1]+1]), len(nomiss), len(accmiss))
    pivot = pd.DataFrame({'Accession': list(missing_dict.keys()), 'NumMissing': nomiss, 'AccMissing': accmiss})
    with open(os.path.join(o, 'stats.tsv'), 'w') as handle:
        pivot.to_csv(handle, header=True, index=False, sep='\t')

if __name__ == '__main__':
    compare()

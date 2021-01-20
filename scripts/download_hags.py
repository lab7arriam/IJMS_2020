#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import chain

## process accession list

hag_list = pd.read_csv(os.path.join(*[os.getcwd(), 'hag_full.txt']), header=None, sep='\t')

## introduce yourself to the NCBI server
Entrez.email, Entrez.tool = 'yu.malovichko@arriam.ru', 'MyCustomScript'

## write a function to fetch NCBI Protein references

def fetchProteins(query:str, geneNames:str, baseFrom:str='nucleotide', baseTo:str='protein'):

    name_list = geneNames.replace('and ', '').split(', ')
    result_list = []

    with Entrez.efetch(db=baseFrom, id=query, rettype='gb', retmode='txt') as handle:
        result = SeqIO.read(handle, 'gb')
        for feature in result.features:
           if feature.type == 'source':
               id = feature.qualifiers.get('strain')[0]
           if feature.type == 'CDS':
                if feature.qualifiers.get('gene') and feature.qualifiers.get('gene')[0] in name_list:
                    sequence = Seq(feature.qualifiers.get('translation')[0])
                    name = feature.qualifiers.get('gene')[0]
                    entry = SeqRecord(seq=sequence, id=id, name=name, description=name)
                    result_list.append(entry)

    return result_list


all_results = list()
for i in range(len(hag_list)):
    result = fetchProteins(hag_list.iloc[i,0], hag_list.iloc[i,1])
    all_results.append(result)

all_results = list(chain.from_iterable(all_results))


name_set = set(map(lambda x: x.name, all_results))

if 'results' not in os.listdir():
    subprocess.run(['mkdir', 'results'])

for name in name_set:
    gene_list = filter(lambda x: x.name == name, all_results)
    with open(os.path.join(*[os.getcwd(), 'results', f'{name}.fasta']), 'w') as handle:
        SeqIO.write(gene_list, handle, 'fasta')

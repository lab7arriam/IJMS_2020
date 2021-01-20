#! /usr/bin/env python3

import click
import requests
import pandas as pd
import urllib
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def fetchSequences(email, input, output):

    Entrez.email, Entrez.tool = email, 'MyCustomScript'

    output_list = []

    dige_df = pd.read_csv(input, header=0, sep='\t')
    accessions = set(dige_df.iloc[:,0].to_list())
    print(dige_df)
    for acc in accessions:
        try:
            with Entrez.efetch(db='protein', id=acc, rettype='fasta', retmode='txt') as handle:
                result = SeqIO.read(handle, 'fasta')
                result.id = acc
                output_list.append(result)
        except urllib.error.HTTPError:
            url = f'https://www.uniprot.org/uniprot/{acc}.fasta'
            request = requests.get(url, allow_redirects=True)
            request = request.content.decode('utf-8')
            request_meta = request.split('\n')[0]
#            id = request_meta.split(' ')[0].lstrip('>')
            description = ' '.join(request_meta.split(' ')[1:])
            seq = Seq(''.join(request.split('\n')[1:]))
            print(acc)
            record = SeqRecord(id=acc, seq=seq, description=description)
            output_list.append(record)

    with open(output, 'w') as handle:
        SeqIO.write(output_list, handle, 'fasta')

@click.command()
@click.option('--email', '-e', type=str)
@click.argument('input', type=click.File())
@click.argument('output', type=click.Path(exists=False))

def fetch(email:str, input:str, output:str):

    fetchSequences(email, input, output)

if __name__ == '__main__':
    fetch()

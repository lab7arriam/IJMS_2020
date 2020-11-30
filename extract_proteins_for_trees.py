#/usr/bin/python3.7

import argparse
import os
import csv
from collections import defaultdict
from Bio import SeqIO

def update_defaultdict(defdict, key, feature):
    """
    Updates defaultdict with special features 
    If value is empty, creates new list, else appends list
    """
    if key not in defdict:
        defdict[key]=[feature]
    else:
       defdict[key].append(feature)
    return(defdict)


def parse_proteome(gene_dict, proteome: str, extracted_dict):
    """
      Extracts records enlisted in gene_dict from proteome fasta_file
      Takes dictionary with extracted entries and extends it 
    """
    for record in SeqIO.parse(proteome,"fasta"):
        for gene in gene_dict:
            if record.id in gene_dict[gene]:
                record.id=filename.replace('.proteome.faa','')
                update_defaultdict(extracted_dict, gene, record)
    return(extracted_dict)



def create_gene_dict(parse_list):
    """
        Creates gene dict from clustered proteins in roary output
        Key - gene symbol, value - list of protein records
    """
    gene_dict=defaultdict(list)
    for string in parse_list:
        gene_dict[string[0]]=string[1:]
    return(gene_dict)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract proteins from pangenome')
    parser.add_argument('-pd', '--prt_dir', dest='prot_dir', help='path to the proteomes',
                        type=str)
    parser.add_argument('-cls', '--clust_st', dest='clust_stat', help='protein clusters',type=argparse.FileType('r'))
    args = parser.parse_args()

    parse_list=[row.replace('\n','').split() for row in args.clust_stat]
    gene_dict=create_gene_dict(parse_list)
    extracted_dict=defaultdict(list)

    for filename in os.listdir(args.prot_dir):
        if '.faa' in filename: 
            parse_proteome(gene_dict, filename, extracted_dict)
    
    for gene in extracted_dict:
        SeqIO.write(extracted_dict[gene],gene+'.fasta',"fasta")

        


    
    
    
   

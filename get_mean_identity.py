#/usr/bin/python3.7

import argparse
import os
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio import pairwise2 as pw2
from itertools import combinations
from statistics import mean
import multiprocessing
from annotate_bt_assemblies import write_csv

def find_sequence_identity(seq_tuple):
    """
      Performes global pairwise alignmet of two sequences
      Returns blast-like identity between sequences
    """

    global_align = pw2.align.globalxx(seq_tuple[0], seq_tuple[1])
    seq_length = min(len(seq_tuple[0]), len(seq_tuple[1]))
    matches = global_align[0][2]
    percent_match = (matches / seq_length) * 100

    return(percent_match)

def parse_fasta_dir(fasta_dir, ret_type='table', threads=5):
    """
      Performes global pairwise alignmet of two sequences
      Returns blast-like identity between sequences
    """

    seq_dict=defaultdict(list)
    for filename in os.listdir(fasta_dir):
        if 'msa' in filename:
            for record in SeqIO.parse(os.path.join(fasta_dir, filename),"fasta"):
                 seq_dict[filename.replace('.msa','')].append(str(record.seq).replace('-',''))

    pool = multiprocessing.Pool(threads)

    if ret_type=='table':
        ret_list=[['tree','identity']]
        for key in seq_dict:
           print('counting', key)
           id_list=pool.map(find_sequence_identity, list(combinations(seq_dict[key][0:2], 2)))
           ret_list.append([key, mean(id_list)])
        return(ret_list)

    elif ret_type=='matrix':
        pairs=list()
        heatmap=[['']+list(seq_dict.keys())]

        for key1 in seq_dict:
            key1_dict={}
            print(key1)
            for key2 in seq_dict:
                if key1==key2:
                    if len(seq_dict[key1])==1:
                        key1_dict[key1]=100
                    else:
                        key1_dict[key1]=mean(pool.map(find_sequence_identity, list(combinations(seq_dict[key1], 2))))
                else:
                    key1_dict[key2]=mean(pool.map(find_sequence_identity, list(combinations([seq_dict[key1]+seq_dict[key2]][0], 2))))
            heatmap.append([key1]+[key1_dict[key] for key in seq_dict.keys() if key in key1_dict])
        return(heatmap)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Counts mean identity of the fasta files')
    parser.add_argument('-md', '--m_dir', dest='msa_dir', help='path to the MSA files',
                        type=str)
    parser.add_argument('-t', '--th', dest='threads', help='the number of threads',
                        type=str, default=5)
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()
    
    write_csv(parse_fasta_dir(args.msa_dir),args.out_file)
    

        


    
    
    
   

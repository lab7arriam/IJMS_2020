#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from collections import defaultdict
from annotate_bt_assemblies import read_csv_to_list, create_dict_from_list, write_csv


def parse_fasta_files(analyzed_dir,flagellin_dict, *extr_acc : str):
    """
        Searches for fasta files in the specified dir and extracts fasta sequences
        Returns a list with metadata
    """
    genomes_dict=defaultdict(dict)
    count_dict=defaultdict(dict)
    for filename in os.listdir(analyzed_dir):
        for record in SeqIO.parse(os.path.join(analyzed_dir, filename),"fasta"):
            if filename.split(':')[0] not in genomes_dict[record.id]:
                genomes_dict[record.id][filename.split(':')[0]]=record
                count_dict[record.id][filename.split(':')[0]]=1
            else:
                count_dict[record.id][filename.split(':')[0]]+=1
                genomes_dict[record.id][filename.split(':')[0]+'.'+str(count_dict[record.id][filename.split(':')[0]])]=record
                flagellin_dict[filename.split(':')[0]+'.'+str(count_dict[record.id][filename.split(':')[0]])]=flagellin_dict[filename.split(':')[0]]
    filtered_seqs=[]
    ret_list=[['assembly', 'gene_list', 'length_list', 'max_gene','max_len', 'max_num']]
    if extr_acc:
        extr_dict=defaultdict(dict)
        for acc in extr_acc:
           extr_dict[acc]=[]
    for assembly in genomes_dict:
         if extr_dict:
             if sum([acc in genomes_dict[assembly].keys() for acc in extr_dict.keys() ])==2:
                 for acc in  extr_dict.keys():
                   extr_dict[acc].append(genomes_dict[assembly][acc])
         #creating indicies list corresponding to the adundance
         sort_inds=sorted([ int(flagellin_dict[key][0]) for key in genomes_dict[assembly].keys()],reverse=True)
         for i in sort_inds:
             for key in genomes_dict[assembly].keys():
                 if int(flagellin_dict[key][0])==i:
                     #choosing the maximum key
                     max_key=key
             if len(genomes_dict[assembly][max_key].seq)>200:
                 #checking the sequence length
                 break        
         filtered_seqs.append(genomes_dict[assembly][max_key])
         ret_list.append([assembly,'; '.join([ key for key in genomes_dict[assembly].keys()]), 
                        '; '.join([str(len(genomes_dict[assembly][key].seq)) for key in genomes_dict[assembly].keys()]),
                         max_key,len(genomes_dict[assembly][max_key].seq), i])
    #SeqIO.write(filtered_seqs,'flagellin_refs.fasta',"fasta")
    #for acc in extr_dict:
    #    SeqIO.write(extr_dict[acc],acc+'.fasta',"fasta")
    return(ret_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts flagelin sequences for the Bt assemblies')
    parser.add_argument('-f', '--f_dir', dest='fast_dir', help='the directory with the extracted fasta files',
                        type=str)
    parser.add_argument('-c', '--csv_t', dest='csv_tab', help='the roary gene presence absence result',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()
  
    flagellin_dict=create_dict_from_list(read_csv_to_list(args.csv_tab, headless=True,delim=','),0,3)
    write_csv(parse_fasta_files(args.fast_dir,flagellin_dict,'hag_1','BK713_14900'),args.out_file)

    



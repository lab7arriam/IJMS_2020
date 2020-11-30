#/usr/bin/python3.7

import argparse
import os
from Bio import SeqIO
from collections import defaultdict
from annotate_bt_assemblies import read_csv_to_list, write_csv

def summarize_diamond_dat(diamond_res, acc_list):
    """
      Parses diamond-generated tsv-table and chooses the hit with the lowes p-value
      Returns the list with gene presence-absence according to the order accessions in the gene_list
    """
    diamond_res_dict=defaultdict(dict)
    for row in read_csv_to_list(diamond_res, headless=False,delim='\t'):
        diamond_res_dict[row[1]][row[0]]=[row[2],row[3],row[10]]
    ret_list=[diamond_res.split('/')[-1].split('.')[0]]
    for acc in acc_list:
        if acc not in list(diamond_res_dict.keys()):
            ret_list.append('0')
        else:
            passed_dict=defaultdict(list)
            for gene in diamond_res_dict[acc]:
                if float(diamond_res_dict[acc][gene][2]) <= min([float(diamond_res_dict[acc][gene][2]) for gene in diamond_res_dict[acc]]) and float(diamond_res_dict[acc][gene][0])>67:
                    passed_dict[gene]=diamond_res_dict[acc][gene]
            if len(passed_dict)<1:
                ret_list.append('0')
            if len(passed_dict)>=1:
                second_pass=defaultdict(list)
                for gene in passed_dict:
                    if float(passed_dict[gene][0]) >= max([float(passed_dict[gene][0]) for gene in passed_dict]):
                        second_pass[gene]=passed_dict[gene]
                for gene in second_pass:
                    if float(second_pass[gene][1]) >= max([float(second_pass[gene][1]) for gene in second_pass]):
                        ret_list.append(gene)
                        break
    
    return(ret_list)
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='summarizes gene presence/absence data based on diamond blastp results')
    parser.add_argument('-t', '--t_dir', dest='tsv_dir', help='the directory with the diamond results',
                        type=str)
    parser.add_argument('-a', '--a_tab', dest='acc_tab', help='a table with accession numbers of the analysed proteins',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()
    accessions=[item[0] for item in read_csv_to_list(args.acc_tab, headless=True,delim=',')]
    #print(len(accessions))
    screened_results=[['assembly']+accessions]
    for diamond_tab in os.listdir(args.tsv_dir):
        screened_results.append(summarize_diamond_dat(os.path.join(os.path.realpath(args.tsv_dir),diamond_tab),accessions))
    write_csv(screened_results,args.out_file)



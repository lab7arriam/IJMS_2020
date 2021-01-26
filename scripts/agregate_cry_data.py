#/usr/bin/python3.7

import argparse
import os
from collections import defaultdict
from annotate_bt_assemblies import read_csv_to_list, create_dict_from_list, write_csv



def parse_fasta_files(annot_dict : defaultdict(list), cry_dir : str):
    """
        Parses CryProcessor results in the directory, looks at the diamond_matches.txt files
        Returns a list with agregated information about cry toxins for each assembly
    """

    cry_prots_dict=defaultdict(list)
    for assembly in annot_dict:
        cry_prots_dict[assembly.split('.')[0]]=[annot_dict[assembly][0],annot_dict[assembly][2]]
    cry_test_dict=defaultdict(list)

    for dirname in os.listdir(cry_dir):
        assembly=dirname
        dirname=os.path.join(cry_dir,dirname, 'logs')

        filt_list=list()
        for filename in os.listdir(dirname):
            if 'diamond_matches' in filename:
                parse_list=read_csv_to_list(os.path.join(dirname, filename),headless=False)
                for parse_hit in parse_list:
                    if float(parse_hit[2])>99:
                        filt_list.append(parse_hit[1])
                    else:
                        filt_list.append(parse_hit[1] +'_'+ parse_hit[2])
        cry_prots_dict[assembly].append('; '.join(filt_list))

    ret_list=[cry_prots_dict[assembly] for assembly in cry_prots_dict.keys()]
    ret_list.insert(0, ['assembly', 'serovar', 'cry_list'])
    return(ret_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Agregates the results of CryProcessor on analysed bt assemblies')
    parser.add_argument('-c', '--c_dir', dest='cry_dir', help='the directory with the CryProcessor results',
                        type=str)
    parser.add_argument('-a', '--an_tab', dest='annot_tab', help='the table with assemblies\' annotation',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()
  
    annot_dict=create_dict_from_list(read_csv_to_list(args.annot_tab, headless=True,delim='\t'),0)
    
    write_csv(parse_fasta_files(annot_dict, args.cry_dir),args.out_file)

    



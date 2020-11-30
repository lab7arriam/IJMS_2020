#/usr/bin/python3.7

import argparse
from collections import defaultdict
from annotate_bt_assemblies import read_csv_to_list, create_dict_from_list, write_csv


def dict_from_txt(raw_list, clean_flag=False):
    """
        Reads parsed raw table with flagellin annotations, which contains strings with assemblies names and gff strings
        Uses the assembly name as a key in the output dictionary
    """
    out_flag_dict=defaultdict(list)
    for an_list in raw_list:
        if len(an_list)==1:
            key=an_list[0]
            out_flag_dict[key]=[[],[]]
        else:
            try:
                if not clean_flag:
                    out_flag_dict[key][0].append(str(int(an_list[8].split('length.')[1])+1 ))
                    out_flag_dict[key][1].append(an_list[8].split('product=')[1].split(';')[0])
                else:
                    if 'cap' not in an_list[8] :
                        out_flag_dict[key][0].append(str(int(an_list[8].split('length.')[1])+1 ))
                        out_flag_dict[key][1].append(an_list[8].split('product=')[1].split(';')[0])
            except:
                if not clean_flag:
                    out_flag_dict[key][0].append('pseudo')
                    out_flag_dict[key][1].append('pseudo')
                else:
                    pass
    return(out_flag_dict)

def comapre_dicts(ref_dict,raw_dict):
    """
        Compares values from the dictionaries with flagellin data
        Returns a list with rows to be saved as a csv-file
    """
    ret_list=[['assembly', 'roary_genes','roary_lens', 'roary_num', 'def_genes', 'def_lens', 'def_num']]
    for key in ref_dict:
        ret_list.append([key, ref_dict[key][0],ref_dict[key][1], len(ref_dict[key][1].split('; ')),'; '.join(raw_dict[key][1]), '; '.join(raw_dict[key][0]), len(raw_dict[key][0])])
    return(ret_list)
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compares the number of flagellin genes for annotated gffs')
    parser.add_argument('-r', '--r_tab', dest='roar_tab', help='the roary based flagellin data',
                        type=str)
    parser.add_argument('-d', '--def_t', dest='def_tab', help='the raw text file from the default assemblies',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()
  
    flagellin_dict=create_dict_from_list(read_csv_to_list(args.roar_tab, headless=True,delim='\t'),0,1,2)
    raw_flag_dict=dict_from_txt(read_csv_to_list(args.def_tab, headless=False,delim='\t'),clean_flag=True) 
    write_csv(comapre_dicts(flagellin_dict,raw_flag_dict),args.out_file)
    

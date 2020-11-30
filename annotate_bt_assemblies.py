#/usr/bin/python3.7

import argparse
import csv
from collections import defaultdict


def create_dict_from_list(parse_list, key_ind, *val_inds):
    """
        Creates list dict from parsing list
        Key refers to the asserted index (key_ind), value - by default takes all the indexes (including the key index)
    """
    parse_dict=defaultdict(list)
    for string in parse_list:
        if not val_inds:
            parse_dict[string[key_ind]]=string
        else:
            parse_dict[string[key_ind]]=[string[i] for i in range(len(string)) if i in val_inds]
    return(parse_dict)


def read_csv_to_list(in_file, headless=True, delim='\t'):
    """
        Reads csv file and returns list without header by default 
        If headless argument is false, parses the whole file
    """
    ret_list=list()
    with open(in_file,'r') as csv_file:
        my_reader = csv.reader(csv_file, delimiter=delim) 
        if headless:
            next(my_reader)
        for row in my_reader:
            ret_list.append(list(row))
    return(ret_list)


def write_csv(row_list,out_name,*header_strings : str):
    """
       A universal function for writing lists to csv-files
       If input is list of lists uses writerows function else iteratively writes 
       If strings for header are specified, writes header to the output, saves headless table otherwise
    """
    with open(out_name,'w',newline='') as result_file:
        wr = csv.writer(result_file, delimiter='\t')
        if header_strings:
            wr.writerow([name for name in header_strings])
        if type(row_list[0]) is list:
            wr.writerows(row_list)
        else:
            for row in row_list:
                wr.writerow([row])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract proteins from pangenome')
    parser.add_argument('-ba', '--bac_a', dest='bac_as', help='path to the assembly table',
                        type=str)
    parser.add_argument('-st', '--str_t', dest='str_tbl', help='path to the strains table',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()

    assembly_dict=create_dict_from_list(read_csv_to_list(args.bac_as, headless=False),4,1,2,4)
    strain_dict=create_dict_from_list(read_csv_to_list(args.str_tbl, headless=False),0,1)

    for key in assembly_dict:
        assembly_dict[key][0]=assembly_dict[key][0].split(' serovar ')[1].split()[0]
        assembly_dict[key].append(strain_dict[key][0])
        assembly_dict[key]=[assembly_dict[key][i] for i in [2,1,0,3]]
    write_csv([assembly_dict[key] for key in assembly_dict.keys()],args.out_file,'GenBank_accession','Assembly_level', 'Serovar', 'Strain')
    


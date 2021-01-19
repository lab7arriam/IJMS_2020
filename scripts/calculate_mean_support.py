#/usr/bin/python3.7

import csv
import argparse
from statistics import mean


def get_mean_support_for_tree(tree_file, sup_type='r'):
    """
        Calculates mean support for the tree in newick format
        If type is regular (r, by deafult), integer support numbers are expected, otherwise parses tree according to the expected posision of the supporting values
    """
    tree_str=''
    with open (tree_file, 'r',newline='') as csvfile2:
        my_reader3 = csv.reader(csvfile2, delimiter='\t')
        for row in my_reader3:
            tree_str=row[0]
    if sup_type=='r':
        tree_list=tree_str.replace('(',',').replace(':',',').replace(')',',').split(',')
        parsed_list=[]

        for el in tree_list:
            if '_' not in el and '.' not in el and len(el)>0 and ';' not in el:
                try:
                    parsed_list.append(int(el)) 
                except:
                    pass

        return(mean(parsed_list))
    else:
        parsed_list=[]
        tree_list=tree_str.replace('(','').replace(',',')').split(')')
        for el in tree_list:
            if '_' not in el and len(el)>0 and ';' not in el and el.split(':')[0]!='0':
                parsed_list.append(float(el.split(':')[0]))
        return(mean(parsed_list))
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculates mean supporting value for the tree in the newick format')
    parser.add_argument('-t', '--tr_f', dest='tree_file', help='path to the analyzed tree file',
                        type=str)
    args = parser.parse_args()
    try:
        print(get_mean_support_for_tree(args.tree_file,sup_type='i'))
    except:
        print(get_mean_support_for_tree(args.tree_file,sup_type='r'))

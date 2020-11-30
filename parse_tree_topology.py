#/usr/bin/python3.7

import argparse
from collections import defaultdict
from annotate_bt_assemblies import read_csv_to_list, write_csv
from ete3 import Tree
from typing import List


def create_assemblies_dict(parse_list):
    """
        Creates list dict from parsing list
        Key refers to the serovar, value - the list of assemblies
    """
    parse_dict=defaultdict(list)
    for string in parse_list:
        parse_dict[string[2]].append(string[0])
    return(parse_dict)

def rename_nodes(tree):
    for node in tree:
        node.name=node.name.split('.')[0]+'.'+node.name.split('.')[1].split('_')[0]
    return(tree)

def gca(tree: Tree, nodes: List[str]) -> Tree:
    no = [tree & n for n in nodes]
    
    if len(no) == 1:
        return no[0]
    return tree.get_common_ancestor(no)

def get_ancestors(strains_dist, tree: Tree, full_flag=True):
    ret_list=[]
    for key in strains_dist:
        if len(strains_dist[key])>1:
            if full_flag:
                ret_list.append([key, len(strains_dist[key]), len(gca(tree,strains_dist[key]))])
            else:
                ret_list.append(len(gca(tree,strains_dist[key])))
    return(ret_list)


def merge_lists(list1,list2):
    comb_list=[]
    for i in range(len(list1)):
        comb_list.append(list(list1[i])+[list2[i]])
    return(comb_list)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Checks the consistency of serovars clustering in the phylogenetic tree')
    parser.add_argument('-d', '--d_tab', dest='dat_tab', help='the tsv-table describing strains and assemblies',
                        type=str)
    parser.add_argument('-t', '--tree', dest='tree_file', help='trees to analysze', type=argparse.FileType('r'),
                        action='append')
    parser.add_argument('-o', '--out', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()
    strains_dist=create_assemblies_dict(read_csv_to_list(args.dat_tab, headless=True,delim='\t'))
    trees = []
    for tf in args.tree_file:
        t = Tree(tf.name)
        trees.append(rename_nodes(t))
    rows_list=[]
    rows_list=get_ancestors(strains_dist, trees[0])
    for tree in trees[1:]:
        num_list=get_ancestors(strains_dist, tree,full_flag=False)
        rows_list=merge_lists(rows_list,num_list)
    write_csv(rows_list,args.out_file,'serovar','def_len','core', 'bin', 'mash','mini','gyr', 'nprb','guaB','sucC','dnaK','mmsA','yjlD',
'rph','groEL','atpD','rocA','fusA','phbB','tuf','inhA','calY','flagellin')
    #write_csv(rows_list,args.out_file,'serovar','def_len','prot_untr', 'prot_tr', 'nucl_untr','nucl_tr','PCR_untr','PCR_tr')

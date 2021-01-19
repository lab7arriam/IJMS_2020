#/usr/bin/python3.7

import argparse
import os
import subprocess
from annotate_bt_assemblies import read_csv_to_list,create_dict_from_list
from Bio import SeqIO
from collections import defaultdict
import itertools
from statistics import mean
from sklearn.feature_extraction import DictVectorizer
import numpy as np
import sys
sys.setrecursionlimit(25000)


def run_minimap(analyzed_dir, reference, query, out_dir, out_index):

    """
        Runs minimap2 aligner with the reference and the query specified
    """
   
    cmd_dist = subprocess.call('cd {0}; mkdir {1}; minimap2 -x asm5  {2} {3} --secondary=no --paf-no-hit  -o {4} '.format(os.path.realpath(analyzed_dir), out_dir, reference, query, os.path.join(out_dir, out_index)),shell=True)


def process_fastas(assembly_dict, analyzed_dir,out_dir, run_flag=True):

    """
        Prepares input for the minipam aligner comparisions 
 	The most complete reference is chosen as a reference; if the assemly levels are equal, the assembly with the lowest number of contigs is prefered
        Returns the dictionary enshrining the assemblies' data
    """

    comparision_dict={'Complete Genome': 5,'Chromosome': 4,'Scaffold': 3,'Contig': 2}
    analysed_dict=defaultdict(list)

    for filename in os.listdir(analyzed_dir):
        if 'fna' in filename:
           key=filename.split('.')[0]+'.'+filename.split('.')[1].split('_')[0]
           analysed_dict[key]=[len(list(SeqIO.parse(os.path.join(analyzed_dir, filename),"fasta"))),comparision_dict[assembly_dict[key][0]],filename, sum([len(record.seq) for record in list(SeqIO.parse(os.path.join(analyzed_dir, filename),"fasta"))])]

    combinations = list(itertools.combinations(list(analysed_dict.keys()),2))
    refs=[]

    for pair in combinations:

        if analysed_dict[pair[0]][1] > analysed_dict[pair[1]][1]:
            refs.append((pair[0],pair[1]))

        elif analysed_dict[pair[0]][1] < analysed_dict[pair[1]][1]:
            refs.append((pair[1],pair[0]))

        else:
            if analysed_dict[pair[0]][0] < analysed_dict[pair[1]][0]:
                refs.append((pair[0],pair[1]))

            elif analysed_dict[pair[0]][0] < analysed_dict[pair[1]][0]:
                refs.append((pair[1],pair[0]))

            else:
                refs.append((pair[0],pair[1]))

    if run_flag:
        for pair in refs:
            run_minimap(analyzed_dir, analysed_dict[pair[0]][2],analysed_dict[pair[1]][2],out_dir,pair[0]+':'+pair[1]+'.paf')

    return(analysed_dict)

def getOverlap(a, b,check=True):

    """
        Returns the overlap between the intervals
        By default returns interception if check flag is set to False returns interval in the following format: ((end, begin),len)
    """

    if check:
        return max(0, min(a[0][0], b[0][0]) - max(a[0][1], b[0][1]))
    else:
        if a[0][0]==b[0][0] and a[0][1]==b[0][1]:

            coeff=max(float(a[2]),float(b[2]))
            match_len=int(coeff*int((max(a[0][0], b[0][0]) - min(a[0][1], b[0][1]))))

        else:
            len1,len2,sum_len, over_len=a[0][0]-a[0][1], b[0][0]-b[0][1], max(a[0][0], b[0][0])-min(a[0][1], b[0][1]),max(0, min(a[0][0], b[0][0]) - max(a[0][1], b[0][1]))
 
            coeff=a[2]*(len1-over_len)/(sum_len)+b[2]*(len2-over_len)/(sum_len)+over_len/sum_len*((a[2]+b[2])/2)
            match_len=int(sum_len*coeff)

            if sum_len*coeff < max(a[1],b[1]):

               match_len=(int((a[1]+b[1]+max(a[1],b[1]))/2))
               coeff=match_len/sum_len
       
        return(((max(a[0][0], b[0][0]),min(a[0][1], b[0][1])),match_len,coeff))



def check_overlaps(interval_list):

    """
        Checks if any of the intervals presented in the list oveplap
        Returns True if an overlap exists and False otherwise
    """

    for pair1 in interval_list:
        for pair2 in interval_list:
            if pair1!=pair2:
                if int(getOverlap(pair1,pair2))>0:
                    return(True)
                    break
    return(False)


def reduce_interval_list(interval_list):

    """
        Recursively reduces the list with the intervals
        If an overlap exists returns the union of the intervals
    """
    reducing_pairs=list(set(interval_list))
    ret_set=set()
    #re-write code excluding nested cycles
    for pair1 in reducing_pairs:
        flag=0
        for pair2 in reducing_pairs:
            if pair1!=pair2:
       
                if int(getOverlap(pair1,pair2))>0:
                    flag=1
                    ret_set.add(getOverlap(pair1,pair2,check=False))
                    break
        if flag==0:
            ret_set.add(pair1)

    if check_overlaps(list(ret_set)):
        return(reduce_interval_list(list(ret_set)))
    return(list(ret_set))


def overlap_intervals(interval_dict):

    """
        Iterates over the intervals dictionary and returns pooled intervals
        Launches reduce_interval_list()
    """

    ret_dict=defaultdict(list)

     #use set to discard duplicate intervals
    for cont in interval_dict:
        interval_dict[cont]=list(set(interval_dict[cont]))

    for cont in interval_dict:
     
        if len(interval_dict[cont])>1:
            ret_dict[cont]=reduce_interval_list(interval_dict[cont])

        # if is the contig has only one interval we do not perform anything
        else:
            ret_dict[cont]=interval_dict[cont]

    return(ret_dict)



def count_mean_id(paf_list,count_type='sets',allow_zeroes=True):

    """
        Parses the PAF-table's row list and calculates the mean identity between the genomes
    """

    if count_type=='contigs':
        #counts the mean blast-like identity for each match
        ids=[]

        for row in paf_list:
            if row[4]=='*':
                #appends non-mapped contigs if
 the allow_zeroes flag is specified 
                if allow_zeroes:
                    ids.append(0)
                else:
                    pass
            else:
                ids.append(int(row[9])/int(row[10]))

        return(mean(ids))

    elif count_type=='sets':
        #counts the double length of the matched unionited intervals

        sets_dist=defaultdict(list)
        for row in paf_list:
            if row[4]!='*':
                try:
                    sets_dist[row[5]].append(((int(row[8]),int(row[7])),int(row[9]),int(row[9])/int(row[10])))
                except:
                    pass
        overlapped=overlap_intervals(sets_dist)
 
        return(sum([sum([interval[1] for interval in overlapped[contig]]) for contig in overlapped])*2)
        

def agregate_minimap_results(analyzed_dict,analyzed_dir,out_dir):

    """
        Iterates over the PAF results of the genomes minimap2-generated genome alignments
        Returns a numpy array with the results
    """
    resulting_dict=dict()

    for filename in os.listdir(os.path.join(analyzed_dir,out_dir)):
        
        ref_set,query_set=set(),set()
        paf_file=read_csv_to_list(os.path.join(os.path.join(analyzed_dir,out_dir), filename), headless=False)
        
        #find the contig lists for the query and the target
        for row in paf_file:
            ref_set.add(row[5])
            query_set.add(row[0])
        
        #if contigs are absent, append the row specifying unmapped hits
        if len(query_set) < analyzed_dict[filename.split(':')[1].replace('.paf','')][0]:
            paf_file.extend([['-']*4+['*']+['-']*5 for j in range(int(analyzed_dict[filename.split(':')[1].replace('.paf','')][0])-len(query_set))])
        if len(ref_set) < analyzed_dict[filename.split(':')[0]][0]:
            paf_file.extend([['-']*4+['*']+['-']*5 for j in range(int(analyzed_dict[filename.split(':')[0]][0])-len(ref_set))])
       
        #sum(alignments)*2/sum(genomes)
        print(filename.replace('.paf',''))

        match_len = count_mean_id(paf_file)
        identity = match_len/(int(analyzed_dict[filename.split(':')[1].replace('.paf','')][3])+int(analyzed_dict[filename.split(':')[0]][3]))

        if identity>1:
            identity=match_len/2/max(int(analyzed_dict[filename.split(':')[1].replace('.paf','')][3]),int(analyzed_dict[filename.split(':')[0]][3]))
            if identity>1: 
                identity = count_mean_id(paf_file,count_type='contigs')
                if identity>1: 
                    identity = 1-min(int(analyzed_dict[filename.split(':')[1].replace('.paf','')][3]),int(analyzed_dict[filename.split(':')[0]][3]))/max(int(analyzed_dict[filename.split(':')[1].replace('.paf','')][3]),int(analyzed_dict[filename.split(':')[0]][3]))

        resulting_dict[filename.replace('.paf','')]=identity

    print(resulting_dict)
          

    i=0
    data_dict = [{} for i in range(len(analyzed_dict))]
    
    for key1 in sorted(analyzed_dict.keys()):
        for key2 in sorted(analyzed_dict.keys()):

            if key1+':'+key2 in resulting_dict:
                data_dict[i][analyzed_dict[key2][2]]=resulting_dict[key1+':'+key2]
            elif key2+':'+key1 in resulting_dict:
                data_dict[i][analyzed_dict[key2][2]]=resulting_dict[key2+':'+key1]
            elif key1==key2:
                data_dict[i][analyzed_dict[key2][2]]=1
        i+=1

    dictvectorizer = DictVectorizer(sparse=False)
    features = dictvectorizer.fit_transform(data_dict)

    feature_names = dictvectorizer.get_feature_names()
    
    rows = np.array(feature_names, dtype='<U80')[:, np.newaxis]
    feature_names.insert(0,'')
    

    cols = np.array(feature_names, dtype='<U80')
    features=np.vstack((cols,np.hstack((rows, features))))

    return(features)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs genome comparision with minimap2 and creates a comparision matrix of the results')
    parser.add_argument('-f', '--fa_dir', dest='fast_dir', help='path to the fasta files',
                        type=str, default=None)
    parser.add_argument('-c', '--c_tab', dest='csv_tab', help='table with the properties of genomes',
                        type=str)
    parser.add_argument('-od', '--o_dir', dest='out_dir', help='the name of the output directory',
                        type=str, default='Alignments')
    parser.add_argument('-o', '--o_file', dest='out_file', help='the name of the output file',
                        type=str)
    args = parser.parse_args()

    assembly_dict = create_dict_from_list(read_csv_to_list(args.csv_tab, headless=False),0,1)
    result_dict=process_fastas(assembly_dict, args.fast_dir,args.out_dir,run_flag=True)
    dist_array= agregate_minimap_results(result_dict, args.fast_dir,args.out_dir)

    with open(args.out_file, 'w') as f:
        np.savetxt(f, dist_array, delimiter='\t', fmt='%s')

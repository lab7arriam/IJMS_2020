#! /usr/bin/env python3

import click
import os
import pandas as pd
from Bio import SeqIO

@click.command()
@click.option('-t', type=str)
@click.option('-i', type=str)
@click.option('-c', type=str)
@click.option('-o', type=str)

def pivot(t, i, c, o):

    with open(t, 'r') as handle:
        template = pd.read_csv(handle, header=0, sep='\t')
        template.assembly = pd.Series(map(lambda x: x.split('.')[0], template.assembly))

    with open(c, 'r') as handle2:
        all_clusters = pd.read_csv(handle2, header=0, sep=',', keep_default_na=False)
        assemblies_names = [x for x in all_clusters.columns.values if x.find('GCA') > -1]
        all_clusters_dict = dict()
        for assembly in assemblies_names:
            local_list = []
            for j in range(0, all_clusters.shape[0]):
                if all_clusters[assembly][j] != '':
                    local_list.append(all_clusters[all_clusters.columns.values[0]][j])
            all_clusters_dict[assembly.split('.')[0]] = ';'.join(local_list)

    cluster_dict = dict()
    files = os.listdir(i)
    for file in files:
        prefix = '_'.join(file.split('_')[0:2])
        cluster_set = set()
        with open(os.path.join(i, file), 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                description = record.description.split(' ')
                cluster = list(filter(lambda x: x.find('cluster') > -1, description))[0]
                cluster = cluster.split('=')[1]
                cluster_set.add(cluster)
        cluster_string = ';'.join(cluster_set)
        cluster_dict[prefix] = cluster_string

    sorter = list(cluster_dict.keys())

    template.assembly = template.assembly.astype('category')
    template.assembly.cat.set_categories(sorter, inplace=True, ordered=True)
    template.sort_values(['assembly'], inplace=True)

    all_clusters_dict = {k: all_clusters_dict[k] for k in sorter}

    template.insert(1, 'Initial_clusters', list(all_clusters_dict.values()), True)
    template.insert(3, 'PCR_products', list(cluster_dict.values()), True)

    template.to_csv(o, header=True, index=False, sep='\t')

if __name__ == '__main__':
    pivot()

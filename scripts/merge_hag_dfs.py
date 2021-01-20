#! /usr/bin/env python3

import click
import os
import pandas as pd

my_df = pd.read_csv(os.path.join('pivots', 'hag_search_stats_1e-10.tsv'), header=0, sep='\t')
not_my_df = pd.read_csv(os.path.join('pivots', 'flagellin_dat_lengths.csv'), header=0, sep='\t')
lengths = pd.read_csv(os.path.join('pivots', 'hag_default_length_1e-10.tsv'),header=0, sep='\t')

sorter = my_df.Accession.to_list()
#print(sorter)

not_my_df.assembly = list(map(lambda x: x.split('.')[0], not_my_df.assembly))
not_my_df.length_list = list(map(lambda x: x.replace(' ', ''), not_my_df.length_list))
not_my_df.assembly = not_my_df.assembly.astype('category')
not_my_df.assembly.cat.set_categories(sorter, inplace=True, ordered=True)
not_my_df.sort_values(['assembly'], inplace=True)
#print(not_my_df.assembly)

lengths.Accession = lengths.Accession.astype('category')
lengths.Accession.cat.set_categories(sorter, inplace=True, ordered=True)
lengths.sort_values(['Accession'], inplace=True)
#print(lengths.Default_length.to_list())

my_df.insert(2, 'Default_lengths', lengths.Default_length.to_list(), True)
#print(not_my_df.length_list)
my_df.insert(4, 'Annot_lengths', not_my_df.length_list.to_list(), True)
my_df.insert(4, 'Annot', list(map(lambda x: x.count(';') + 1, not_my_df.length_list.to_list())), True)

my_df.to_csv(os.path.join('pivots', 'hag_full_stat_upd_flag_1e-10.tsv'), header=True, index=False, sep='\t')

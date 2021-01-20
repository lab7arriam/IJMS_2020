#! /usr/bin/env python3

import click
import os
import pandas as pd

@click.command()
@click.option('--input', '-i', help='Path to the directory containing tables to merge', type=click.Path(exists=True))
@click.option('--output', '-o', help='Path to the output file', type=click.Path(exists=False))

def lengthPivot(input, output):

    sorter = []
    df_dict = dict()
    files = os.listdir(input)
    for file in files:
        prefix = file.split('_')[2]
        with open(os.path.join(input, file), 'r') as handle:
            df = pd.read_csv(handle, sep='\t', header=None)
            df[df.columns.values[0]] = pd.Series(map(lambda x: '_'.join(x.split('_')[0:2]),  df[df.columns.values[0]]))
            if len(sorter) == 0:
                sorter = sorted(df[df.columns.values[0]].to_list())
            df_dict[prefix] = df

    out_df = pd.DataFrame({'Accession': sorter})
    for key in sorted(df_dict.keys()):
        df = df_dict[key]
        tosort = df.columns.values[0]
        df[tosort] = df[tosort].astype('category')
        df[tosort].cat.set_categories(sorter, inplace=True)
        df.sort_values([tosort], inplace=True)

        pos = out_df.shape[1]
        out_df.insert(pos, key, df[df.columns.values[2]].to_list(), True)

    with open(output, 'w') as handle:
        out_df.to_csv(handle, sep='\t', header=True, index=False)

if __name__ == '__main__':
    lengthPivot()

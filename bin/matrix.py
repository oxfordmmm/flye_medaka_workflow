#!/usr/bin/env python3
import pandas as pd
from argparse import ArgumentParser


def run(opts):
    dfs = []
    for i in opts.input:
        df = pd.read_csv(i)
        dfs.append(df)

    df = pd.concat(dfs)

    #meta = pd.read_csv(opts.meta)
    #meta_dict = dict(zip(meta['Barcode'], meta['Sample']))
    #print(meta_dict)
    #df['Sample name 1'] = df['Sample name 1'].map(meta_dict)
    #df['Sample name 2'] = df['Sample name 2'].map(meta_dict)
    df=df[df['type']=='SNPs']
    df2=df.pivot_table(index='Sample name 1', values='Ref', columns='Sample name 2', aggfunc='first')
    df2.to_csv(opts.output)

if __name__ == '__main__':
    parser = ArgumentParser(description='Parse dna diff and add meta data')
    parser.add_argument('-i', '--input', required=True, nargs='+',
                        help='input csvs files')
    parser.add_argument('-m', '--meta', required=False,
                        help='meta data')
    parser.add_argument('-o', '--output', required=True,
                        help='output tsv file')


    opts, unknown_args = parser.parse_known_args()

    run(opts)
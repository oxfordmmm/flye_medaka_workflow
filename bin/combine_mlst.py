#!/usr/bin/env python3
import pandas as pd
from argparse import ArgumentParser


def run(opts):
    dfs = []
    for i in opts.input:
        if opts.scheme == '':
            df = pd.read_csv(i, sep='\t', names=['seqid', 'scheme', 'ST', 'gene 1', 'gene 2', 'gene 3', 'gene 4', 'gene 5', 'gene 6', 'gene 7', 'gene 8', 'gene 9', 'gene 10'])
        else:
            df = pd.read_csv(i, sep='\t')
        dfs.append(df)

    df = pd.concat(dfs)

    #meta = pd.read_csv(opts.meta)
    #meta['seq_name'] = meta['Barcode'].astype(str) + '.fasta'
    #meta_dict = dict(zip(meta['seq_name'], meta['Sample']))
    #print(meta_dict)
    #df['Sample name'] = df['seqid'].map(meta_dict)
    
    if opts.scheme == '':
        cols=['Sample name', 'seqid', 'scheme', 'ST', 'gene 1', 'gene 2', 'gene 3', 'gene 4', 'gene 5', 'gene 6', 'gene 7', 'gene 8', 'gene 9', 'gene 10']
        df=df[cols]
    df.sort_values(by='FILE', inplace=True)
    df.to_csv(opts.output, index=False)

if __name__ == '__main__':
    parser = ArgumentParser(description='Parse dna diff and add meta data')
    parser.add_argument('-i', '--input', required=True, nargs='+',
                        help='input csvs files')
    parser.add_argument('-m', '--meta', required=True,
                        help='meta data')
    parser.add_argument('-o', '--output', required=True,
                        help='output tsv file')
    parser.add_argument('-s', '--scheme', required=False, default='',
                        help='MLST scheme used')


    opts, unknown_args = parser.parse_known_args()

    run(opts)
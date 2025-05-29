#!/usr/bin/env python3
from Bio import SeqIO
import sys

def filter_contigs(input_file, output_file, min_length):
    for record in SeqIO.parse(input_file, 'fasta'):
        if len(record.seq) >= min_length:
            contig_name=record.id
            with open(f'{contig_name}_{output_file}', 'w') as out:
                record.id = sys.argv[4] + '_' + record.id
                SeqIO.write(record, out, 'fasta')

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_length = int(sys.argv[3])
    filter_contigs(input_file, output_file, min_length)
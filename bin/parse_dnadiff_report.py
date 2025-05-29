#!/usr/bin/env python3
import pandas as pd
import re
import six
from collections import defaultdict
from collections import namedtuple
from argparse import ArgumentParser

Property = namedtuple('Property', 'ref query')
PropertyWithPerc = namedtuple('PropertyWithPerc', 'ref ref_perc query query_perc')
def _parse_dnadiff_into_sections(report_file):
    """Parse dnadiff output lines into sections."""
    report_fh = open(report_file, 'r')
    section = "NO_SECTION"
    sections = defaultdict(list)
    for line in report_fh:
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('/') or line.startswith('NUCMER') or line.startswith('[REF]'):
            continue
        if line.startswith('['):
            section = line
            section = section.replace('[', '')
            section = section.replace(']', '')
        else:
            sections[section].append(line)
    return sections

def _parse_percent_field(field):
    """Parse dnadiff field with percent value."""
    tmp = field.split('(')
    perc = tmp[1].replace(')', '')
    perc = perc.replace('%', '')
    return float(tmp[0]), float(perc)

def _parse_simple_section(lines):
    """Parse a simple dnadiff report section."""
    results = {}
    for line in lines:
        tmp = re.split("\s+", line)
        if '%' not in tmp[1] and '%' not in tmp[2]:
            results[tmp[0]] = Property(float(tmp[1]), float(tmp[2]))
        else:
            ref_prop, ref_prop_perc = _parse_percent_field(tmp[1])
            query_prop, query_prop_perc = _parse_percent_field(tmp[2])
            results[tmp[0]] = {'ref prop': ref_prop, 
                    'ref prop %': ref_prop_perc, 
                    'query prop': query_prop,
                    'query prop %': query_prop_perc}
            #PropertyWithPerc(ref_prop, ref_prop_perc, query_prop, query_prop_perc)
    return results


def _parse_complex_section(lines):
    """Parse a complex dnadiff report section."""
    section = "NO_SECTION"
    sections = defaultdict(list)
    results = defaultdict(dict)
    # Parse alignment section into subsections:
    for line in lines:
        if len(line) == 0:
            continue
        # FIXME: Very specific to current dnadiff output:
        if line.startswith('1-to-1') or line.startswith('M-to-M') or re.match("Total(S|G|I)", line):
            tmp = re.split("\s+", line)
            section = tmp[0]
            results[section]['Number'] = {'Ref':tmp[1], 'Query':tmp[2]}# Property(float(tmp[1]), float(tmp[2]))
        else:
            sections[section].append(line)

    # Parse subsections and update results dictionary:
    for section, lines in six.iteritems(sections):
        parsed = _parse_simple_section(lines)
        for name, prop in six.iteritems(parsed):
            results[section][name] = prop
    return results


def parse_file(f):
    sections=_parse_dnadiff_into_sections(f)
    results_snps = _parse_complex_section(sections['SNPs'])

    results_snps['TotalIndels']['Number']['type']='INDELs'
    results_snps['TotalSNPs']['Number']['type']='SNPs'

    return results_snps['TotalIndels']['Number'],results_snps['TotalSNPs']['Number']

def run(opts):
    l=[]
    indels,snps=parse_file(opts.input)
    l.append(indels)
    l.append(snps)
    df=pd.DataFrame(l)
    df['Sample name 1']=opts.sample_name1
    df['Sample name 2']=opts.sample_name2
    #print(df)
    df.to_csv(opts.output, index=False)


if __name__ == '__main__':
    parser = ArgumentParser(description='Parse dna diff and add meta data')
    parser.add_argument('-i', '--input', required=True,
                        help='input TSV file')
    parser.add_argument('-s1', '--sample_name1', required=True,
                        help='sample_name 1')
    parser.add_argument('-s2', '--sample_name2', required=True,
                        help='sample_name 2')
    parser.add_argument('-o', '--output', required=True,
                        help='output tsv file')


    opts, unknown_args = parser.parse_known_args()

    run(opts)

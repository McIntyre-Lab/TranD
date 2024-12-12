#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:43:51 2024

@author: k.bankole
"""

import argparse
import pandas as pd


def getOptions():

    parser = argparse.ArgumentParser(
        description='Add genome as prefix and ERP as suffix for jxnHash in input GTF file.')

    # Define the arguments
    parser.add_argument(
        '-i',
        '--inGTF',
        required=True,
        dest="inGTF",
        help='Path to the GTF file.',
    )

    parser.add_argument(
        '-e',
        '--erpFile',
        required=True,
        dest="erpFile",
        help='Path to the ERP file.'
    )

    parser.add_argument(
        '-g',
        '--genome',
        required=True,
        dest="genome",
        help='Genome identifier (e.g., dmel6).'
    )

    parser.add_argument(
        '-o',
        '--outfile',
        required=True,
        dest="outfile",
        help='Output file path for the new gtf.'
    )

    # Parse arguments and return them
    return parser.parse_args()


def main():
    inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    erpFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er_vs_fiveSpecies_2_dmel6_ujc_infoERP.csv"
    genome = 'dmel6'
    outfile = ""

    inGTF = args.inGTF
    erpFile = args.erpFile
    genome = args.genome
    outfile = args.outfile

    columns = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']

    data = pd.read_csv(inGTF, sep='\t', comment='#',
                       header=None, names=columns,
                       low_memory=False)

    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
    data['geneID'] = data['attribute'].str.extract(r'gene_id "([^"]+)"')
    data['transcriptID'] = data['attribute'].str.extract(
        r'transcript_id "([^"]+)"')

    erpDfr = pd.read_csv(erpFile, low_memory=False)[
        ['jxnHash', 'ERP']].set_index('jxnHash')

    erpDct = erpDfr['ERP'].to_dict()

    data['transcriptID'] = data['transcriptID'].apply(
        lambda x: genome + "_" + x + "_" + erpDct[x])

    data['attributes'] = data.apply(
        lambda row: f'transcript_id "{row["transcriptID"]}"; gene_id "{row["geneID"]}";', axis=1)
    data.drop(columns=['geneID', 'transcriptID'], inplace=True)

    gtf_columns = ['seqname', 'source', 'feature', 'start',
                   'end', 'score', 'strand', 'frame', 'attributes']
    data = data[gtf_columns]

    # Output the DataFrame to GTF format
    data.to_csv(outfile, sep='\t', header=False, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

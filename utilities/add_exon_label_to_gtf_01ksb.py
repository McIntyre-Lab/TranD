#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io
import csv

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Add a label for the number of each exon feature per transcript")

    # Input data
    parser.add_argument(
        "-g",
        "--gtf",
        dest="gtfFile",
        required=True,
        help="Path to GTF"
    )
    
    parser.add_argument(
        "-l",
        "--label",
        dest="label",
        default="exon_",
        required=True,
        help="Name of label (label_1, label_2, ...). Default: exon"
    )

    # Output data
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Path and name of outfile"
    )

    args = parser.parse_args()
    return args

def main():
    
    # gtfFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"
    # label = "ER"
    # outfile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/labelled_fiveSpecies_2_dmel6_ujc_er.gtf"
    
    gtfFile = args.gtfFile
    label = args.label
    outfile = args.outfile
        
    dfr = trand.io.read_exon_data_from_file(gtfFile)
    
    dfr['exonLabel'] = dfr.groupby('transcript_id').cumcount() + 1
    dfr['exonLabel'] = label + '_' + dfr['exonLabel'].astype(str)
        
    dfr.loc[:, 'source'] = "TranD"
    dfr.loc[:, 'feature'] = "exon"
    dfr.loc[:, 'score'] = "."
    dfr.loc[:, 'frame'] = "."
    # dfr.loc[:, 'attribute'] = dfr.apply(lambda row: 'transcript_id "{}"; gene_id "{}"; {}_number "{}";'.format(row['transcript_id'],row['gene_id'],label,row['exonLabel']), axis=1)
    dfr.loc[:, 'attribute'] = dfr.apply(lambda row: 'transcript_id "{}_{}"; gene_id "{}";'.format(row['transcript_id'], row['exonLabel'], row['gene_id']), axis=1)
                                         
    output_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                           'frame', 'attribute']
    
    dfr = dfr.reindex(columns=output_column_names)
    dfr.to_csv(outfile, sep="\t", mode='a', index=False, header=False,
                doublequote=False, quoting=csv.QUOTE_NONE)




if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

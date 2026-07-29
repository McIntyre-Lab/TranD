#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io
import csv


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Convert a GTF into gene_model format. "
        "A gene model GTF contains one transcript per gene, and "
        "each exon is considered an exon region "
        "and no exon region should be overlapping another within the same "
        "gene.")

    # Input data
    parser.add_argument(
        "-g",
        "--gtf",
        dest="gtfFile",
        required=True,
        help="Path to GTF"
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
    # outfile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/labelled_fiveSpecies_2_dmel6_ujc_er.gtf"

    # gtfFile = "C:/Users/knife/Desktop/fiveSpecies_2_dmel6_ujc_er.gtf"
    # gtfFile = "C:/Users/knife/Desktop/dmel-all-r6.50.gtf"
    # outfile = "C:/Users/knife/Downloads/test.gtf"

    gtfFile = args.gtfFile
    outfile = args.outfile

    dfr = trand.io.read_exon_data_from_file(gtfFile)
    # dfr = read_exon_data_from_file(gtfFile)

    # Check every gene has one transcript
    if not dfr.groupby('gene_id')['transcript_id'].nunique().eq(1).all():
        raise Exception("There are genes with more than one transcript. "
                        "This GTF cannot be a gene model.")
    else:
        print("Every gene has one transcript. Continuing...")

    # Check for overlapping exons

    dfr = dfr.sort_values(['seqname', 'gene_id', 'start'], ignore_index=True)
    records = dfr.to_dict('records')

    prevRow = None
    for row in records:
        if prevRow == None:
            prevRow = row
        else:
            prevSeqname = prevRow['seqname']
            seqname = row['seqname']

            prevGene = prevRow['gene_id']
            gene = row['gene_id']

            if prevSeqname != seqname:
                # Dont compare between chromosomes
                prevRow = row
            elif prevGene != gene:
                # Dont compare between genes
                prevRow = row
            else:
                prevStart = prevRow['start']
                prevEnd = prevRow['end']

                start = row['start']
                end = row['end']

                if max(prevStart, start) < min(prevEnd, end):
                    print(prevRow, row)
                    raise Exception("The above rows indicate overlapping exons. "
                                    " within a gene. This GTF cannot be a gene "
                                    " model in its current state."
                                    )
                else:
                    prevRow = row
    print("There are no overlapping exons within a gene. This GTF can be a gene model.")

    # Create a column with ER number to include on GTF (for fun)
    dfr['ER'] = dfr.groupby(['seqname', 'gene_id'])[
        ['start', 'end']].cumcount() + 1
    dfr['ER'] = dfr['gene_id'] + '_ER' + dfr['ER'].astype(str)

    print()

    print("Setting both gene_id and transcript_id to be the gene_id...")
    dfr['transcript_id'] = dfr['gene_id']

    print("Converting to TranD gene_model format...")
    dfr.loc[:, 'source'] = "TranD"
    dfr.loc[:, 'feature'] = "exon"
    dfr.loc[:, 'score'] = "."
    dfr.loc[:, 'frame'] = "."
    dfr.loc[:, 'attribute'] = dfr.apply(lambda row: 'transcript_id "{}"; gene_id "{}"; ER "{}"'.format(
        row['transcript_id'], row['gene_id'], row['ER']), axis=1)

    output_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                           'frame', 'attribute']

    dfr = dfr.reindex(columns=output_column_names)
    dfr.to_csv(outfile, sep="\t", index=False, header=False,
               doublequote=False, quoting=csv.QUOTE_NONE)

    print("Complete!")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

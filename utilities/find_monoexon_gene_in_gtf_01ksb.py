#!/usr/bin/env python

import argparse
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Input a GTF. Outputs a list of genes in that GTF with "
                    "only monoexon transcripts (single exon genes).")

    # Input data
    parser.add_argument(
        "-g",
        "--gtf",
        dest="gtfFile",
        required=True,
        help="Input GTF"
    )

    # Output data
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Output file. Directory must already exist."
    )

    args = parser.parse_args()
    return args


def main():

    gtfFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    outfile = "/nfshome/k.bankole/Desktop/test_folder/list_fiveSpecies_2_dmel6_single_exon_gene.csv"
    gtfFile = args.gtfFile
    outfile = args.outfile

    gtfDfr = trand.io.read_exon_data_from_file(gtfFile)

    grpDfr = gtfDfr.groupby(['gene_id', 'transcript_id']).size().reset_index()
    grpDfr.columns = ['geneID', 'transcriptID', 'numExon']
    grpGrpDfr = grpDfr.groupby(['geneID']).agg({'numExon': max}).reset_index()

    singleExonGeneDfr = grpGrpDfr[grpGrpDfr['numExon'] == 1]['geneID']
    singleExonGeneDfr.to_csv(outfile, index=False, header=None)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

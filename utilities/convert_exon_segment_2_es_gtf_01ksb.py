#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io
import os


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert event analysis file from TranD 1 GTF gene output "
                                                 "(event_analysis_er.csv) to a GTF where each transcript is "
                                                 "represents all of the exon regions of a gene.")

    # Input data
    parser.add_argument("-e",
                        "--event-analysis",
                        dest="eaFile",
                        required=True,
                        help="Path to event analysis file")

    # Output data
    parser.add_argument("-o",
                        "--output-gtf",
                        dest="outGTF",
                        required=True,
                        help="Path and filename for output GTF")

    args = parser.parse_args()
    return args


def main():

    # eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_EA2GTF/event_analysis_er.csv"
    # gtfOutPath = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_EA2GTF/test_gtf.gtf"

    esFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_1gtf_geneMode_fiveSpecies_sexDetSubset/fiveSpecies_2_dmel6_ujc_sexDetSubset_event_analysis_ef.csv"
    gtfOutPath = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/TEST_fiveSpecies_2_dmel6_ujc_sexDetSubset_es.gtf"

    # eaFile = args.eaFile
    # gtfOutPath = args.outGTF

    esDf = pd.read_csv(esFile, low_memory=False)

    esDf = esDf.rename(columns={'ef_chr': 'seqname',
                                'ef_strand': 'strand',
                                'ef_start': 'start',
                                'ef_end': 'end'})

    esDf['transcript_id'] = esDf['gene_id']

    esDf = esDf[['seqname', 'start', 'end',
                 'strand', 'transcript_id', 'gene_id']]

    outExonDf = esDf.sort_values(
        by=['seqname', 'transcript_id', 'start']).reset_index(drop=True)
    outExonDf['start'] = outExonDf['start'].apply(lambda x: 1 if x == 0 else x)

    if os.path.isfile(gtfOutPath):
        os.remove(gtfOutPath)

    trand.io.write_gtf(data=outExonDf, out_fhs={
                       "gtf": gtfOutPath}, fh_name="gtf")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Output the number of exon segments in a gene that "
                    "are below or equal to an input threshold.")

    # Input data
    parser.add_argument("-i",
                        "--in-GTF",
                        dest="inGTF",
                        required=True,
                        help="Input ES GTF")

    parser.add_argument("-t",
                        "--threshold",
                        dest="threshold",
                        required=True,
                        help="Exon Segment Length Threshold")

    # Output data
    # parser.add_argument("-", "--", dest="", required=True, help="")

    args = parser.parse_args()
    return args


def main():

    inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dser1_ujc_es.gtf"
    threshold = 1

    inGTF = args.inGTF
    threshold = args.threshold

    threshold = int(threshold)
    gtfDf = trand.io.read_exon_data_from_file(inGTF)

    gtfDf['ES'] = 'ES' + (gtfDf.groupby('gene_id').cumcount() + 1).astype(str)
    gtfDf['length'] = gtfDf['end'] - gtfDf['start']

    gtfDf['flag_lngthBlwThresh'] = (gtfDf['length'] <= threshold).astype(int)

    gnWThreshDf = gtfDf.groupby(
        'gene_id')['flag_lngthBlwThresh'].sum().reset_index()
    gnWThreshDf = gnWThreshDf.rename(
        columns={'flag_lngthBlwThresh': 'numSgmntBlwThrshld'})

    print("There are {} segments of length {} (or shorter).".format(
        gnWThreshDf['numSgmntBlwThrshld'].sum(), threshold))

    numGnWThreshDf = gnWThreshDf[gnWThreshDf['numSgmntBlwThrshld']
                                 >= 1]['gene_id'].nunique()

    pctGnWThreshDf = numGnWThreshDf / len(gnWThreshDf)

    print("There are {} genes with at least one segment below the threshold. This accounts for {:0.2%} of the total genes in the annotation.".format(
        numGnWThreshDf, pctGnWThreshDf))

    # TODO: ADD OUTPUT IF DESIRED


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

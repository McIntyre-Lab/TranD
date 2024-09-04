#!/usr/bin/env python

import argparse

import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Count num reads per jxnHash of a catted "
                                     "UJC xscript link file.")

    # Input data
    parser.add_argument(
        "-l",
        "--linkFile",
        dest="linkFile",
        required=True,
        help="Catted link file (no multiGene UJC)"
    )

    # Output data
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output directory. Must already exist."
    )

    parser.add_argument(
        "-x",
        "--prefix",
        dest="prefix",
        required=True,
        help="Output prefix. Required."
    )

    args = parser.parse_args()
    return args


def main():

    linkFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_xscript_link_noMultiGene.csv"

    outdir = "whaaat"
    prefix = "huhhhhh"

    linkFile = args.linkFile
    outdir = args.outdir
    prefix = args.prefix

    nmgLinkDf = pd.read_csv(linkFile, low_memory=False)

    nmgLinkDf = nmgLinkDf.rename({'source': 'sampleID'})

    countDf = nmgLinkDf.groupby(['source', 'jxnHash', 'geneID']).count()[
        'transcriptID'].reset_index()

    countDf = countDf.rename(
        columns={'transcriptID': 'numRead'})
    countDf = countDf[['source', 'geneID', 'jxnHash', 'numRead']]

    readPerSampleID = countDf.groupby('source').count()[
        'numRead'].reset_index()
    readPerGene = countDf.groupby('geneID').count()['numRead'].reset_index()

    outPrefix = outdir + f"/{prefix}"

    countDf.to_csv(
        f"{outPrefix}_ujc_count.csv", index=False)

    readPerSampleID.to_csv(
        f"{outPrefix}_read_per_sampleID.csv", index=False)

    readPerGene.to_csv(
        f"{outPrefix}_read_per_gene.csv", index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

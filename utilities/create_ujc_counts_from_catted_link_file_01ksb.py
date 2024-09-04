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

    print("Link file read complete...")

    nmgLinkDf = nmgLinkDf.rename(columns={'source': 'sampleID'})

    countDf = nmgLinkDf.groupby(['sampleID', 'jxnHash', 'geneID']).count()[
        'transcriptID'].reset_index()

    countDf = countDf.rename(
        columns={'transcriptID': 'numRead'})
    countDf = countDf[['sampleID', 'geneID', 'jxnHash', 'numRead']]

    print("Initial count complete...")

    readPerSampleID = countDf.groupby(
        'sampleID')['numRead'].sum().reset_index()

    print("SampleID count complete...")

    readPerGene = countDf.groupby('geneID')['numRead'].sum().reset_index()

    print("Gene count complete...")

    readPerJxnhash = countDf.groupby(
        'jxnHash')['numRead'].sum().reset_index()

    print("Jxnhash count complete...")

    outPrefix = outdir + f"/{prefix}"

    countDf.to_csv(
        f"{outPrefix}_ujc_count.csv", index=False)

    readPerSampleID.to_csv(
        f"{outPrefix}_read_per_sampleID.csv", index=False)

    readPerGene.to_csv(
        f"{outPrefix}_read_per_gene.csv", index=False)

    readPerJxnhash.to_csv(
        f"{outPrefix}_read_per_jxnHash.csv", index=False)

    print("Output complete!")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

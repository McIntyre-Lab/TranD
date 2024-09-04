#!/usr/bin/env python

import argparse

import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Add flagMultiGene to catted xscript link "
                                     "file using ujc dsc file. Split xscript "
                                     "link file into multiGene and nonMultiGene "
                                     "and output.")

    # Input data
    parser.add_argument(
        "-d",
        "--dscFile",
        dest="dscFile",
        required=True,
        help="Catted ujc dscrptn file that has gone through "
        "find_multiGene_jxnHash_in_catted_ujc_dsc_file_01ksb.py")

    parser.add_argument(
        "-l",
        "--linkFile",
        dest="linkFile",
        required=True,
        help="Catted UJC link file")

    # Output data
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output directory. Must already exist.")

    parser.add_argument(
        "-x",
        "--prefix",
        dest="prefix",
        required=True,
        help="Output prefix. Required.")

    args = parser.parse_args()
    return args


def main():

    dscFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_dscrptn.csv"
    linkFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_xscript_link_simpleCat.csv"

    outdir = "yipee"
    prefix = "wahoo"

    dscFile = args.dscFile
    linkFile = args.linkFile
    outdir = args.outdir
    prefix = args.prefix

    dscDf = pd.read_csv(dscFile, low_memory=False)
    linkDf = pd.read_csv(linkFile, low_memory=False)

    print("Num unique jxnHash:", dscDf['jxnHash'].nunique())
    print("Num multiGene jxnHash:",
          dscDf[dscDf['flagMultiGene'] == 1]['jxnHash'].nunique())
    print("Num nonMultiGene jxnHash:",
          dscDf[dscDf['flagMultiGene'] == 0]['jxnHash'].nunique())
    print("Num total read:", len(linkDf))

    dscForMergeDf = dscDf[['jxnHash', 'flagMultiGene']].drop_duplicates()

    linkWFlagDf = pd.merge(linkDf, dscForMergeDf, how='outer',
                           on='jxnHash', indicator='merge_check')

    if (linkWFlagDf['merge_check'] != 'both').any():
        raise Exception("An error occurred with the merge. There are jxnHashes "
                        "in the dscrptn file that are not in the link file or  "
                        "vice versa.")
    else:
        linkWFlagDf.drop('merge_check', axis=1, inplace=True)

    nonMultiGeneDf = linkWFlagDf[linkWFlagDf['flagMultiGene'] == 0].copy()
    multiGeneDf = linkWFlagDf[linkWFlagDf['flagMultiGene'] == 1].copy()

    print("Num multiGene jxnHash post-merge:",
          multiGeneDf['jxnHash'].nunique())
    print("Num multiGene reads:", multiGeneDf['transcriptID'].nunique())

    pctMultiGene = len(multiGeneDf) / len(linkWFlagDf)

    print(f"Percent of all READS that are multiGene: {pctMultiGene:%}")

    outPrefix = outdir + f"/{prefix}"

    nonMultiGeneDf.to_csv(
        f"{outPrefix}_ujc_xscript_link_noMultiGene.csv", index=False)

    multiGeneDf.to_csv(
        f"{outPrefix}_ujc_xscript_link_multiGene.csv", index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

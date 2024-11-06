#!/usr/bin/env python

import argparse
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Add flagKeepIR to fiveSpecies flag file.")

    # Input data
    parser.add_argument(
        "-l",
        "--list-file",
        dest="inLst",
        required=True,
        help="List of genes missing from event analysis.")

    parser.add_argument(
        "-f",
        "--flag-file",
        dest="inFlg",
        required=True,
        help="fiveSpecies flag file")

    # Output data
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        help="Name and path of output file. Directory must already exist.")

    args = parser.parse_args()
    return args


def main():

    inLst = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/list_noHeader_fiveSpecies_2_dmel6_genes_missing_from_event_analysis.csv"
    inFlg = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/flag_fiveSpecies_2_dmel6_ujc.csv"

    outfile = "/nfshome/k.bankole/Desktop"

    inLst = args.inLst
    inFlg = args.inFlg
    outfile = args.outfile

    gnLstDf = pd.read_csv(inLst, header=None, low_memory=False)
    flgDf = pd.read_csv(inFlg, low_memory=False)

    gnLstDf.columns = ['geneID']

    mergeDf = pd.merge(gnLstDf, flgDf, how='outer', on=[
                       'geneID'], indicator='merge_check')

    if (mergeDf['merge_check'] == 'left_only').any():
        raise Exception("An error ocurred. Please check your gene list.")
    else:
        mergeDf['flagKeepIR'] = mergeDf['merge_check'].apply(
            lambda x: 1 if x == "both" else 0)
        mergeDf.drop('merge_check', axis=1, inplace=True)

    mergeDf.to_csv(outfile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

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

    gnLstDfr = pd.read_csv(inLst, header=None, low_memory=False)
    flgDfr = pd.read_csv(inFlg, low_memory=False)

    gnLstDfr.columns = ['geneID']

    mergeDfr = pd.merge(gnLstDfr, flgDfr, how='outer', on=[
        'geneID'], indicator='merge_check')

    if (mergeDfr['merge_check'] == 'left_only').any():
        raise Exception("An error ocurred. Please check your gene list.")
    else:
        mergeDfr['flagKeepIR'] = mergeDfr['merge_check'].apply(
            lambda x: 1 if x == "both" else 0)
        mergeDfr.drop('merge_check', axis=1, inplace=True)

    moveCol = mergeDfr.pop('geneID')
    mergeDfr.insert(1, 'geneID', moveCol)

    mergeDfr.to_csv(outfile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

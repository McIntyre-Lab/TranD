#!/usr/bin/env python

import argparse
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Input concatenated UJC description file "
                                     "from ujc run split byGene. Outputs new "
                                     "dscrptn file with flagMultiGene and no "
                                     "duplicate jxnHash. ")

    # Input data
    parser.add_argument(
        "-i",
        "--catDscFile",
        dest="catDscFile",
        required=True,
        help="Dsc file after catting split gene UJC output")

    # Output data
    parser.add_argument(
        "-o",
        "--outDscFile",
        dest="outDscFile",
        required=True,
        help="Output file + path")

    args = parser.parse_args()
    return args


def main():

    catDscFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_dscrptn_simpleCat.csv"
    outDscFile = "/nfshome/k.bankole/Desktop/test_folder/catdscasiofgn.csv"

    catDscFile = args.catDscFile
    outDscFile = args.outDscFile

    catDscDf = pd.read_csv(catDscFile, low_memory=False)

    # remove the false flag that is a remnant from catting
    catDscDf = catDscDf.drop('flagMultiGene', axis=1)

    print("Num rows in catted dsc: ", len(catDscDf))
    print("Num unique jxnHash:", catDscDf['jxnHash'].nunique())

    uniqOnJxnHashDf = catDscDf.groupby('jxnHash').agg(set)

    tempCol = uniqOnJxnHashDf['geneID'].apply(
        lambda x: len(x) > 1).astype(int)

    uniqOnJxnHashDf = uniqOnJxnHashDf.applymap(
        lambda x: '|'.join(str(value) for value in x))

    uniqOnJxnHashDf['flagMultiGene'] = tempCol

    pctMultiGene = len(
        uniqOnJxnHashDf[uniqOnJxnHashDf['flagMultiGene'] == 1]) / len(uniqOnJxnHashDf)

    print("Num multiGene jxnHash:", len(
        uniqOnJxnHashDf[uniqOnJxnHashDf['flagMultiGene'] == 1]))

    print(f"Percent of all jxnHash that are multiGene: {pctMultiGene:%}")

    outDf = uniqOnJxnHashDf.reset_index()

    outDf.to_csv(outDscFile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

#!/usr/bin/env python

import argparse
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Placeholder.")

    # Input data
    parser.add_argument("-i", "--in-cnt", dest="cntFile",
                        required=True, help="")

    # Output data
    parser.add_argument("-o", "--out-list", dest="outFile",
                        required=True, help="")

    args = parser.parse_args()
    return args


def main():

    cntFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_count.csv"
    outFile = ""

    cntFile = args.cntFile
    outFile = args.outFile

    cntDf = pd.read_csv(cntFile, low_memory=False)

    cntDf[['species', 'sex', 'techRep']
          ] = cntDf['sampleID'].str.split('_', expand=True)

    cntGrpDf = cntDf.groupby(['jxnHash']).agg({'numRead': sum}).reset_index()

    jxnHashGtTen = cntGrpDf[cntGrpDf['numRead'] > 10]['jxnHash']

    jxnHashGtTen.to_csv(outFile, header=False, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

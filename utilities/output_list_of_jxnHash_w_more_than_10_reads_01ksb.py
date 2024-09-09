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

    DIR = "kopp_lmm_head_data"
    SPECIES = "dser"
    GENOME = "dser1"

    cntFile = f"/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/{DIR}/{SPECIES}_data_2_{GENOME}_ujc_count.csv"
    outFile = f"/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/readsWMoreThan10/list_{SPECIES}_data_2_{GENOME}_jxnHash_w_gt_10_read.txt"

    # cntFile = args.cntFile
    # outFile = args.outFile

    cntDf = pd.read_csv(cntFile, low_memory=False)

    # cntDf[['species', 'sex', 'techRep']
    #       ] = cntDf['sampleID'].str.split('_', expand=True)

    cntGrpDf = cntDf.groupby(['jxnHash']).agg({'numRead': sum}).reset_index()

    jxnHashGtTen = cntGrpDf[cntGrpDf['numRead'] > 10]['jxnHash']

    jxnHashGtTen.to_csv(outFile, header=False, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

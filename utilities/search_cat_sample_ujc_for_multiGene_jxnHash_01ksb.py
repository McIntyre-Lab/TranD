#!/usr/bin/env python

import argparse
import pandas as pd
import glob


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Input a \"ujc xscript link\" file that has been run by "
                    "gene and catted together. Script outputs a list of "
                    "jxnHash that exist in multiple genes and a list of jxnHash that is 'oneToOne'.")

    # Input data
    parser.add_argument("-i",
                        "--infile",
                        dest="inFile",
                        required=True,
                        help="Path to input file.")

    # Output data
    parser.add_argument("-p",
                        "--prefix",
                        dest="prefix",
                        required=True,
                        help="Output prefix. Required."
                        "Output will look like: 'list_{prefix}_oneToOne_jxnHash.csv'.")

    parser.add_argument("-o",
                        "--outdir",
                        dest="outdir",
                        required=True,
                        help="Output directory. Must already exist.")

    args = parser.parse_args()
    return args


def main():

    inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dsan_data_2_dsan1_ujc_xscript_link.csv"
    prefix = "test"
    outdir = "/nfshome/k.bankole/Desktop/test_folder/"

    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dmel650_2_dmel6_ujc_xscript_link.csv"
    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dsim202_2_dsim2_ujc_xscript_link.csv"
    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dsimWXD_2_dsim2_ujc_xscript_link.csv"
    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dsan11_2_dsan1_ujc_xscript_link.csv"
    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dyak21_2_dyak2_ujc_xscript_link.csv"
    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dser11_2_dser1_ujc_xscript_link.csv"

    inFile = args.inFile
    outdir = args.outdir
    prefix = args.prefix

    inDf = pd.read_csv(inFile, low_memory=False)

    hash2GeneDf = inDf[['jxnHash', 'geneID']]

    print("Num total rows:", len(hash2GeneDf['jxnHash']))

    uniqHashDf = hash2GeneDf.groupby(
        'jxnHash').agg({'geneID': set}).reset_index()

    print("Num unique jxnHash:", len(uniqHashDf))

    uniqHashDf['flagMultiGene'] = uniqHashDf['geneID'].apply(
        lambda x: len(x) > 1).astype(int)

    multiGeneDf = uniqHashDf[uniqHashDf['flagMultiGene'] == 1].copy()
    oneToOneDf = uniqHashDf[uniqHashDf['flagMultiGene'] == 0].copy()

    pctMultiGene = len(multiGeneDf) / len(uniqHashDf)

    print("Num multiGene jxnHash:", len(multiGeneDf))
    print(f"Percent of all jxnHash that are multiGene: {pctMultiGene:%}")

    outPrefix = outdir + f"/list_{prefix}"

    oneToOneDf['geneID'] = oneToOneDf['geneID'].apply(lambda x: next(iter(x)))
    oneToOneDf = oneToOneDf.drop('flagMultiGene', axis=1)
    oneToOneDf.to_csv(outPrefix + "_oneToOne_jxnHash.csv", index=False)

    multiGeneDf = multiGeneDf.explode('geneID')
    multiGeneDf = multiGeneDf.drop('flagMultiGene', axis=1)
    multiGeneDf.to_csv(outPrefix + "_multiGene_jxnHash.csv", index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

#!/usr/bin/env python

import argparse
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Input a \"ujc xscript link\" file that has been run by "
                    "gene and catted together. Script outputs a list of "
                    "jxnHash that exist in multiple genes.")

    # Input data
    parser.add_argument("-i",
                        "--infile",
                        dest="inFile",
                        required=True,
                        help="Path to input file.")

    # Output data
    parser.add_argument("-o",
                        "--outfile",
                        dest="outFile",
                        required=True,
                        help="Full path for output file.")

    args = parser.parse_args()
    return args


def main():

    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_species_gtf_combined_genes/dsan_data_2_dsan1_ujc_xscript_link.csv"

    # outFile = "/nfshome/k.bankole/Desktop/test_folder/multigene.txt"

    inFile = args.inFile
    outFile = args.outFile

    inDf = pd.read_csv(inFile, low_memory=False)

    hash2GeneDf = inDf[['jxnHash', 'geneID']]

    print("Num total rows:", len(hash2GeneDf['jxnHash']))

    uniqHashDf = hash2GeneDf.groupby(
        'jxnHash').agg({'geneID': set}).reset_index()

    print("Num unique jxnHash:", len(uniqHashDf))

    uniqHashDf['multiGene'] = uniqHashDf['geneID'].apply(
        lambda x: len(x) > 1)

    multiGeneDf = uniqHashDf[uniqHashDf['multiGene']]
    # singleGeneDf = uniqHashDf[~uniqHashDf['multiGene']]

    pctMultiGene = len(multiGeneDf) / len(uniqHashDf)

    print("Num multiGene jxnHash:", len(multiGeneDf))
    print(f"Percent of all jxnHash that are multiGene: {pctMultiGene:%}")

    outDf = multiGeneDf['jxnHash']

    outDf.to_csv(outFile, index=False, header=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

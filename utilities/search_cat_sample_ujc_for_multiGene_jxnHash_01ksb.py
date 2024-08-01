#!/usr/bin/env python

import argparse
import pandas as pd


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # Input data
    parser.add_argument("-", "--", dest="", required=True, help="")

    # Output data
    parser.add_argument("-", "--", dest="", required=True, help="")

    args = parser.parse_args()
    return args


def main():

    inFile = ""

    inDf = pd.read_csv(inFile, low_memory=False)

    hash2GeneDf = inDf[['jxnHash', 'geneID']]

    hash2GeneDf = hash2GeneDf.groupby(
        'jxnHash').agg({'geneID': set}).reset_index()

    hash2GeneDf['multiGene'] = hash2GeneDf['geneID'].apply(
        lambda x: len(x) > 1)

    multiGeneDf = hash2GeneDf[hash2GeneDf['multiGene']]
    singleGeneDf = hash2GeneDf[~hash2GeneDf['multiGene']]


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

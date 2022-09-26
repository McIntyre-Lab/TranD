#!/usr/bin/env python

import argparse
import pandas as pd
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Subset TranD pairwise transcript distance file by list of genes."
        )
    )

    # Input data
    parser.add_argument(
        "-p",
        "--pairwise-distance",
        dest="TD",
        required=True,
        help=(
            "Pairwise transcript distance TranD output CSV file."
        )
    )
    parser.add_argument(
        "-t",
        "--list-type",
        dest="inType",
        required=False,
        choices=['gene_id'],
        default='gene_id',
        help=(
            "Type of feature in inclusion/exclusion list provided "
            "(currently only gene_id"
        )
    )
    parser.add_argument(
        "-i",
        "--include-list",
        dest="inInclude",
        required=False,
        help=(
            "List of gene IDs to include in output from input "
            "pairwise distance TranD file (no header) - cannot be used with exclusion list"
        )
    )
    parser.add_argument(
        "-e",
        "--exclude-list",
        dest="inExclude",
        required=False,
        help=(
            "List of gene IDs to exclude in output from input "
            "pairwise distance TranD file (no header) - cannot be used with inclusion list"
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output CSV file for subset pairwise distance TranD file."
    )

    args = parser.parse_args()
    return args


def main():
    # Get TranD pairwise distance file
    tdDF = pd.read_csv(args.TD, low_memory=False)

    # Get list ID variable - currently will only be gene_id
    listVar = args.inType
    
    # Check that only one list argument was provided
    if args.inInclude is None and args.inExclude is None:
        print("ERROR: No list provided to subset pairwise file on")
        sys.exit()
    elif args.inInclude is not None and args.inExclude is None:
        include = True
        idDF = pd.read_csv(args.inInclude, names=[listVar])
    elif args.inInclude is None and args.inExclude is not None:
        include = False
        idDF = pd.read_csv(args.inExclude, names=[listVar])
    else:
        print("ERROR: Cannot use both inclusion and exclusion list - only provide one or the other")
        sys.exit()

    if include:
        # If list is for inclusion
        # Extract rows of pairwise distance that have a value contained in list of IDs
        filterDF = tdDF[tdDF[listVar].isin(idDF[listVar])]
    else:
        # If list is for exclusion
        # Extract GTF elements that do NOT have a value contained in list of IDs
        filterDF = tdDF[~tdDF[listVar].isin(idDF[listVar])]
        
    # Output filtered file
    filterDF.to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

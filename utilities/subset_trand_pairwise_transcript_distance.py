#!/usr/bin/env python

import argparse
import pandas as pd
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Subset TranD pairwise transcript distance file by list of genes or transcripts."
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
        choices=['gene_id', 'transcript_id'],
        default='gene_id',
        help=(
            "Type of feature in inclusion/exclusion list provided. "
            "Using transcript_id will check for transcript_id values in BOTH "
            "transcript_1 AND transcript_2 columns if 1 list provided. If two "
            "lists provided for transcript_id then first list will be compared "
            "against transcript_1 and second against transcript_2 of pariwise "
            "file. Default: gene_id."
        )
    )
    parser.add_argument(
        "-i",
        "--include-list",
        dest="inInclude",
        action="append",
        required=False,
        help=(
            "List of gene (or transcript) IDs to include in output from "
            "pairwise distance TranD file. List should have no header. "
            "Cannot be used with exclusion list. The -i argument can be used "
            "twice to give two separate lists of transcript_id values. "
            "The first will be compared to transcript_1 and the second to "
            "transcript_2 in the pairwise file. If no filtering is desired for "
            "transcript_1 then provide the word 'all' for the first -i argument "
            "and similarly for transcript_2 use the word 'all' for the second -i "
            "argument."
        )
    )
    parser.add_argument(
        "-e",
        "--exclude-list",
        dest="inExclude",
        action="append",
        required=False,
        help=(
            "List of gene (or transcript) IDs to exclude in output from "
            "pairwise distance TranD file. List should have no header. "
            "Cannot be used with inclusion list. The -i argument can be used "
            "twice to give two separate lists of transcript_id values. "
            "The first will be compared to transcript_1 and the second to "
            "transcript_2 in the pairwise file."
        )
    )
    parser.add_argument(
        "-n1",
        "--name1",
        dest="name1",
        required=False,
        help=(
            "For 2 GTF pairwise distance file, name used for dataset 1 "
            "that is appended to transcript_1 values "
            "(Default: d1 if two lists provided)."
        )
    )
    parser.add_argument(
        "-n2",
        "--name2",
        dest="name2",
        required=False,
        help=(
            "For 2 GTF pairwise distance file, name used for dataset 2 "
            "that is appended to transcript_2 values "
            "(Default: d2 if two lists provided)."
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
        print("ERROR: No list provided to subset pairwise file on.")
        sys.exit()
    elif args.inInclude is not None and args.inExclude is None:
        include = True
        if len(args.inInclude) == 1:
            idDF = pd.read_csv(args.inInclude[0], names=[listVar])
        elif len(args.inInclude) == 2:
            if listVar == "transcript_id":
                if args.inInclude[0] != "all":
                    t1DF = pd.read_csv(args.inInclude[0], names=[listVar])
                    
                    if args.name1 is not None:
                            t1DF[listVar] = t1DF[listVar] + "_" + args.name1
                else:
                    t1DF = tdDF[["transcript_1"]].rename(columns={"transcript_1": listVar})
                if args.inInclude[1] != "all":
                    t2DF = pd.read_csv(args.inInclude[1], names=[listVar])
                    if args.name2 is not None:      
                            t2DF[listVar] = t2DF[listVar] + "_" + args.name2
                else:
                    t2DF = tdDF[["transcript_2"]].drop_duplicates().rename(columns={"transcript_2": listVar})
            else:
                print("ERROR: -t must be transcript_id when more than one list provided.")
                sys.exit()
        else:
            print("ERROR: A maximum of two lists can be provided as input.")
            sys.exit()
    elif args.inInclude is None and args.inExclude is not None:
        include = False
        if len(args.inExclude) == 1:
            idDF = pd.read_csv(args.inExclude[0], names=[listVar])
        elif len(args.inExclude) == 2:
            if listVar == "transcript_id":
                t1DF = pd.read_csv(args.inExclude[0], names=[listVar])
                if args.name1 is not None:
                        t1DF[listVar] = t1DF[listVar] + "_" + args.name1
                t2DF = pd.read_csv(args.inExclude[1], names=[listVar])
                
                if args.name2 is not None:      
                        t2DF[listVar] = t2DF[listVar] + "_" + args.name2
            else:
                print("ERROR: -t must be transcript_id when more than one list provided.")
                sys.exit()
        else:
            print("ERROR: A maximum of two lists can be provided as input.")
            sys.exit()
    else:
        print("ERROR: Cannot use both inclusion and exclusion list - only provide one or the other")
        sys.exit()

    if include:
        # If list is for inclusion
        # Extract rows of pairwise distance that have a value contained in list of IDs
        if listVar == "gene_id":
            filterDF = tdDF[tdDF[listVar].isin(idDF[listVar])]
        else:
            if len(args.inInclude) == 1:
                # A pair is retained if transcript_1 and transcript_2 are BOTH in the list
                filterDF = tdDF[(tdDF["transcript_1"].isin(idDF[listVar])) & (tdDF["transcript_2"].isin(idDF[listVar]))]
            else:
                # A pair is retained if transcript_1 is in list 1 and transcript_2 is in list 2
                filterDF = tdDF[(tdDF["transcript_1"].isin(t1DF[listVar])) & (tdDF["transcript_2"].isin(t2DF[listVar]))]
    else:
        # If list is for exclusion
        # Extract GTF elements that do NOT have a value contained in list of IDs
        if listVar == "gene_id":
            filterDF = tdDF[~tdDF[listVar].isin(idDF[listVar])]
        else:
            if len(args.inInclude) == 1:
                # A pair is retained if transcript_1 AND transcript_2 are NOT in the list
                filterDF = tdDF[(~tdDF["transcript_1"].isin(idDF[listVar])) & (~tdDF["transcript_2"].isin(idDF[listVar]))]
            else:
                # A pair is retained if transcript_1 is NOT in list 1 and transcript_2 is NOT in list 2
                filterDF = tdDF[(~tdDF["transcript_1"].isin(t1DF[listVar])) & (~tdDF["transcript_2"].isin(t2DF[listVar]))]
        
    # Output filtered file
    filterDF.to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    #main()
    
    # Get TranD pairwise distance file
    tdDF = pd.read_csv(args.TD, low_memory=False)

    # Get list ID variable - currently will only be gene_id
    listVar = args.inType
    
    # Check that only one list argument was provided
    if args.inInclude is None and args.inExclude is None:
        print("ERROR: No list provided to subset pairwise file on.")
        sys.exit()
    elif args.inInclude is not None and args.inExclude is None:
        include = True
        if len(args.inInclude) == 1:
            idDF = pd.read_csv(args.inInclude[0], names=[listVar])
        elif len(args.inInclude) == 2:
            if listVar == "transcript_id":
                if args.inInclude[0] != "all":
                    t1DF = pd.read_csv(args.inInclude[0], names=[listVar])
                    
                    if args.name1 is not None:
                            t1DF[listVar] = t1DF[listVar] + "_" + args.name1
                else:
                    t1DF = tdDF[["transcript_1"]].rename(columns={"transcript_1": listVar})
                if args.inInclude[1] != "all":
                    t2DF = pd.read_csv(args.inInclude[1], names=[listVar])
                    if args.name2 is not None:      
                            t2DF[listVar] = t2DF[listVar] + "_" + args.name2
                else:
                    t2DF = tdDF[["transcript_2"]].drop_duplicates().rename(columns={"transcript_2": listVar})
            else:
                print("ERROR: -t must be transcript_id when more than one list provided.")
                sys.exit()
        else:
            print("ERROR: A maximum of two lists can be provided as input.")
            sys.exit()
    elif args.inInclude is None and args.inExclude is not None:
        include = False
        if len(args.inExclude) == 1:
            idDF = pd.read_csv(args.inExclude[0], names=[listVar])
        elif len(args.inExclude) == 2:
            if listVar == "transcript_id":
                t1DF = pd.read_csv(args.inExclude[0], names=[listVar])
                if args.name1 is not None:
                        t1DF[listVar] = t1DF[listVar] + "_" + args.name1
                t2DF = pd.read_csv(args.inExclude[1], names=[listVar])
                
                if args.name2 is not None:      
                        t2DF[listVar] = t2DF[listVar] + "_" + args.name2
            else:
                print("ERROR: -t must be transcript_id when more than one list provided.")
                sys.exit()
        else:
            print("ERROR: A maximum of two lists can be provided as input.")
            sys.exit()
    else:
        print("ERROR: Cannot use both inclusion and exclusion list - only provide one or the other")
        sys.exit()

    if include:
        # If list is for inclusion
        # Extract rows of pairwise distance that have a value contained in list of IDs
        if listVar == "gene_id":
            filterDF = tdDF[tdDF[listVar].isin(idDF[listVar])]
        else:
            if len(args.inInclude) == 1:
                # A pair is retained if transcript_1 and transcript_2 are BOTH in the list
                filterDF = tdDF[(tdDF["transcript_1"].isin(idDF[listVar])) & (tdDF["transcript_2"].isin(idDF[listVar]))]
            else:
                # A pair is retained if transcript_1 is in list 1 and transcript_2 is in list 2
                filterDF = tdDF[(tdDF["transcript_1"].isin(t1DF[listVar])) & (tdDF["transcript_2"].isin(t2DF[listVar]))]
    else:
        # If list is for exclusion
        # Extract GTF elements that do NOT have a value contained in list of IDs
        if listVar == "gene_id":
            filterDF = tdDF[~tdDF[listVar].isin(idDF[listVar])]
        else:
            if len(args.inInclude) == 1:
                # A pair is retained if transcript_1 AND transcript_2 are NOT in the list
                filterDF = tdDF[(~tdDF["transcript_1"].isin(idDF[listVar])) & (~tdDF["transcript_2"].isin(idDF[listVar]))]
            else:
                # A pair is retained if transcript_1 is NOT in list 1 and transcript_2 is NOT in list 2
                filterDF = tdDF[(~tdDF["transcript_1"].isin(t1DF[listVar])) & (~tdDF["transcript_2"].isin(t2DF[listVar]))]
        
    # Output filtered file
    filterDF.to_csv(args.outFile, index=False)
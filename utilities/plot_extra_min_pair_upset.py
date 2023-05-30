#!/usr/bin/env python

import argparse
import trand.plot_functions as PF
import matplotlib.pyplot as plt
import pandas as pd
import sys


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Plot two minimum distance UpSet plots for 1) extras in dataset 1 and 2) extras in dataset 2.")

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="Input minimum pairwise distance file."
    )
    parser.add_argument(
        "-1",
        "--name1",
        dest="inN1",
        default="d1",
        required=False,
        help="Input name of dataset 1 (Default = d1)."
    )
    parser.add_argument(
        "-2",
        "--name2",
        dest="inN2",
        default="d2",
        required=False,
        help="Input name of dataset 2 (Default = d2)."
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output-directory",
        dest="outDir",
        required=True,
        help="Output directory."
    )

    args = parser.parse_args()
    return args

def main():
    # Import input file and check
    td_data= pd.read_csv(args.inFile, low_memory=False)
    name1 = args.inN1
    name2 = args.inN2

    # Check that input is a minimum pairwise distance file with the expected variable names
    if "flag_min_match_"+name1 not in td_data.columns or "flag_min_match_"+name2 not in td_data.columns:
        print("ERROR: Input minimum pairwise distance file does not contain expected columns, check that names provided are correct")
        sys.exit()
    if len(td_data) == 0:
        print("ERROR: Input minimum pairwise distance file is empty.")
        sys.exit()
    checkCols = ["flag_RMP",
                 "flag_alt_exon",
                 "flag_alt_DA",
                 "flag_IR",
                 "flag_5_var",
                 "flag_3_var",
                 "flag_no_ovlp_nt",
                 "num_nt_noOvlp",
                 "prop_nt_noOvlp"]
    if len([c for c in checkCols if c not in td_data.columns])>0:
        print("ERROR: Input minimum pairwise distance file does not contain expected columns")
        print([c for c in checkCols if c not in td_data.columns])
        sys.exit()

    # Plot extras for dataset 1
    PF.plot_min_pair_AS_upset_nt_box(
        td_data[(td_data["flag_min_match_"+name1]==1)&(td_data["flag_RMP"]==0)],
        name1,
        name2,
        "{}/{}.rtf".format(args.outDir, name1 + "_extra_min_pair_AS_upset_nt_box"),
        reciprocal=False
    )
    plt.savefig("{}/{}.png".format(args.outDir, name1 + "_extra_min_pair_AS_upset_nt_box"), dpi=600, format="png")
    plt.clf()


    # # Plot extras for dataset 2
    PF.plot_min_pair_AS_upset_nt_box(
        td_data[(td_data["flag_min_match_"+name2]==1)&(td_data["flag_RMP"]==0)],
        name1,
        name2,
        "{}/{}.rtf".format(args.outDir, name2 + "_extra_min_pair_AS_upset_nt_box"),
        reciprocal=False
    )
    plt.savefig("{}/{}.png".format(args.outDir, name2 + "_extra_min_pair_AS_upset_nt_box"), dpi=600, format="png")
    plt.clf()    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()



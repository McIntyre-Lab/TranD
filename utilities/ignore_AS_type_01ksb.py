# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:55:09 2023

@author: k.bankole
"""

"""
Create a customizable upset plot using existing TranD output

Can remove different types of Alternative Splicing
"""

import plot_functions as PF
import myio
import argparse
import numpy as np
import pandas as pd


# removing all should just show a grey dot w number of pairs
# default is to show all

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=
                                     "Make customizable plots using TranD ouput files "
                                     "with the option to ignore specific types of AS."
                                     "Default plot will show " "all types of AS, 
                                     "ignoring all wil show the number of pairs."
                                     "Input a trand output file (csv) (--indir), ignore options "
                                     "(--ignore-3, --ignore-5 --ignoreAD, --ignoreAE, 
                                     "--ignoreIR, --ignoreNSNT), and output path (--outdir)."
                                     )
    

    # Input data
    parser.add_argument(
        "-i",
        "--indir",
        dest="indir",
        required=True,
        help="Location of csv file containing transcriptome AS analysis")
    
    #Ignore Options
    parser.add_argument(
            "-i3",
            "--ignore-3",
            dest="ignore_3",
            required=False,
            help="Ignore transcript pairs with 3' variation alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-i5",
            "--ignore-5",
            dest="ignore_5",
            required=False,
            help="Ignore transcript pairs with 5' variation alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iAD",
            "--ignore-AD",
            dest="ignore_AD",
            required=False,
            help="Ignore transcript pairs with alternate donor/acceptor alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iAE",
            "--ignore-AE",
            dest="ignore_AE",
            required=False,
            help="Ignore transcript pairs with alternate exon alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iIR",
            "--ignore-IR",
            dest="ignore_IR",
            required=False,
            help="Ignore transcript pairs with intron retention alternative splicing when creating the Upset Plot"
            )
    
    parser.add_argument(
            "-iNSNT",
            "--ignore-NSNT",
            dest="ignore_NSNT",
            required=False,
            help="Ignore transcript pairs with no shared nucleotides when creating the Upset Plot"
            )
    
    # Output data
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help=(
                "Output directory, created if missing. "
                "Default: current directory."
        )
    )
    
    #idea,,, make the prefix the name of the file + whatever was ignored
    # works because the input is already a TranD output file
    # parser.add_argument(
    #     "-x",
    #     "--prefix",
    #     dest="outPrefix",
    #     required=True,
    #     help=(
    #             "Output prefix for plots. "
    #             "Default: no prefix for 1GTF, 'name1_vs_name2' for 2GTF."
    #     )
    # )
    
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )

    args = parser.parse_args()
    return args

# Should this be here or in plot_functions?
def plot_pairAS
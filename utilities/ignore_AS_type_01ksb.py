# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Tue Jan 31 12:55:09 2023

@author: k.bankole
"""

"""
Create a customizable upset plot using existing TranD output

Can remove different types of Alternative Splicing
"""


import trand.plot_functions as PF
import trand.io
import argparse
import numpy as np
import pandas as pd
from sys import exit
import os
import matplotlib.pyplot as plt


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Make customizable plots using TranD ouput files "
                                     "with the option to ignore up to 4 types of AS."
                                     "Default plot will show all types of AS"
                                     "Input a trand output file (csv) (--indir), ignore options (up to 4) "
                                     "(--ignore-3, --ignore-5 --ignoreAD, --ignoreAE, "
                                     "--ignoreIR, --ignoreNSNT), and output path (--outdir)."
                                     )

    # Input data
    parser.add_argument(
        "-i",
        "--indir",
        dest="indir",
        required=True,
        help="Location of csv file containing transcriptome AS analysis")

    # Ignore Options
    parser.add_argument(
        "-i3",
        "--ignore-3",
        dest="ignore_3",
        action="store_true",
        help="Ignore transcript pairs with 3' variation alternative splicing when creating the Upset Plot"
    )

    parser.add_argument(
        "-i5",
        "--ignore-5",
        dest="ignore_5",
        action="store_true",
        help="Ignore transcript pairs with 5' variation alternative splicing when creating the Upset Plot"
    )

    parser.add_argument(
        "-iAD",
        "--ignore-AD",
        dest="ignore_AD",
        action="store_true",
        help="Ignore transcript pairs with alternate donor/acceptor alternative splicing when creating the Upset Plot"
    )

    parser.add_argument(
        "-iAE",
        "--ignore-AE",
        dest="ignore_AE",
        action="store_true",
        help="Ignore transcript pairs with alternate exon alternative splicing when creating the Upset Plot"
    )

    parser.add_argument(
        "-iIR",
        "--ignore-IR",
        dest="ignore_IR",
        action="store_true",
        help="Ignore transcript pairs with intron retention alternative splicing when creating the Upset Plot"
    )

    parser.add_argument(
        "-iNSNT",
        "--ignore-NSNT",
        dest="ignore_NSNT",
        action="store_true",
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

    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite existing output directory and files within.",
    )

    args = parser.parse_args()
    return args



# input Dataframe (csv), Legend Output File Location, Ignore Column Flags
def plotCustomPlot(inDf, legendOut, ignoreDct):
    """
        Function to plot an upset plot based on TranD output with customizable
        AS options.

        Parameters
        ----------
        in_Df : DATAFRAME
                INPUT DATAFRAME TO PLOT.
        legendOut : STRING
                CONTAINS PATH TO PLACE LEGEND OUTPUT.
        ignoreDct : DICTIONARY
                THE DICTIONARY CONTAINING THE AS TYPES AS A STRING AND WHETHER TO IGNORE THEM.

        Returns
        -------
        None.

    """
        
    pairAS = inDf[
        [
            "gene_id",
            "transcript_1",
            "transcript_2",
            "flag_alt_exon",
            "flag_alt_donor_acceptor",
            "flag_IR",
            "flag_5_variation",
            "flag_3_variation",
            "flag_no_shared_nt",
            "num_nt_diff",
            "prop_nt_diff",
        ]
    ].copy()

    # Ensure number of nt different are int and float values
    pairAS["num_nt_diff"] = pairAS["num_nt_diff"].astype(int)
    pairAS["prop_nt_diff"] = pairAS["prop_nt_diff"].astype(float)

    # assures that an exorbitant amt of nucleotides diff are not included in graph
    # bc nsnt will have 100% difference
    pairAS["num_nt_diff"] = np.where(
        pairAS["flag_no_shared_nt"] == 1,
        np.nan,
        pairAS["num_nt_diff"]
    )

    # Make AS flags boolean values and rename
    # converts data frame to only AS t/f and number of NT diff and prop diff
    xcrptFlagCols = [
        "transcript_1",
        "transcript_2",
        "flag_alt_exon",
        "flag_alt_donor_acceptor",
        "flag_IR",
        "flag_5_variation",
        "flag_3_variation",
        "flag_no_shared_nt",
    ]
    pairAS[xcrptFlagCols] = pairAS[xcrptFlagCols].astype(bool)
    pairAS = pairAS.rename(
        columns={
            "flag_alt_exon": "Alt. Exon",
            "flag_alt_donor_acceptor": "Alt. Donor/Acceptor",
            "flag_IR": "Intron Retention",
            "flag_5_variation": "5' Variation",
            "flag_3_variation": "3' Variation",
            "flag_no_shared_nt": "No Shared NT",
            "num_nt_diff": "# NT Different",
            "prop_nt_diff": "Proportion\nNT Different",
        }
    )

    # Create Column Subset based on user arguments
    aSCols = []

    for column_head, ignore_flag in ignoreDct.items():
        if not ignore_flag[0]:
            aSCols.append(column_head)

    if (len(aSCols) > 1):
        PF.plot_upset(
            pairAS.set_index(aSCols),
            "",
            [
                "# NT Different",
                "Proportion\nNT Different",
            ],
        )
        
    legendText = (
        "Number of transcript pairs with the specified types of alternative "
        "splicing indicated by the black dots below the histogram of transcript "
        "pair counts (n = {} transcript pairs, in {} multi-transcript genes). "
        "Box plots represent the number (blue) and proportion (orange) of "
        "nucleotides different between the pair. Transcript pairs with the "
        "alternative splicing methods requested to be ignored are counted "
        "in the first column with no black dots.".format(
            len(pairAS),
            pairAS["gene_id"].nunique()
        )
    )

    with open(legendOut, "w") as outFile:
        startRtf(outFile)
        outFile.write(
            r"\b Figure. Alternative splicing between pairs of transcripts \b0"
            r" \line {}".format(
                legendText
            )
        )
        endRtf(outFile)


# Builds a dictionary based on user input to customize plot headers
def buildIgnoreDct(args):
    """
        Support function to build an ignore dictionary based on user input from argparse.

        Parameters
        ----------
        args : ARGPARSE ARGUMENTS
                THE ARGS FROM ARGPARSE.

        Returns
        -------
        iDct : DICTIONARY
                THE DICTIONARY CONTAINING THE AS TYPES AS A STRING AND WHETHER TO IGNORE THEM.

        """
    iDct = {"3' Variation": [args.ignore_3], "5' Variation": [args.ignore_5], "Alt. Exon": [args.ignore_AE],
              "Alt. Donor/Acceptor": [args.ignore_AD], "Intron Retention": [args.ignore_IR],
              "No Shared NT": [args.ignore_NSNT]}

    return iDct

# uses args to build file name
def buildPrefix(args, ignoreDct):
        """
        Support function to build a prefix for the output file name.

        Parameters
        ----------
        args : ARGPARSE ARGUMENTS
                THE ARGS FROM ARGPARSE.
        ignoreDct : TYPE
                THE DICTIONARY CONTAINING THE AS TYPES AND WHETHER TO IGNORE THEM.

        Returns
        -------
        prefix : STRING
                A PREFIX FOR THE FILE NAME CONTAINING WHICH AS TYPES WERE IGNORED.

        """
        prefix = ""
        ignoreLst = [args.ignore_3, args.ignore_5, args.ignore_AD, args.ignore_AE, args.ignore_IR, args.ignore_NSNT]
        
        if ignoreLst[0]:
                prefix = prefix + "ignore_3_variation_"
        if ignoreLst[1]:
                prefix = prefix + "ignore_5_variation_"
        if ignoreLst[2]:
                prefix = prefix + "ignore_Alt_Donor_Acceptor_"
        if ignoreLst[3]:
                prefix = prefix + "ignore_Alt_Exon_"
        if ignoreLst[4]:
                prefix = prefix + "ignore_Intron_Retention_"
        if ignoreLst[5]:
                prefix = prefix + "ignore_No_Shared_Nucleotides_"
                
        return prefix

#creates legend text
def startRtf(outFile):
    """
        Support function used to output legend text as RTF file.
        
        Parameters
        ----------
        outFile : FILE
                OUTPUT RTF FILE.

        Returns
        -------
        None.

        """
    outFile.write(
        r"{\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\fswiss\fcharset0 Arial;}}"
    )


def getCitation():
    """
        Support function used to output legend text as RTF file.

        Returns
        -------
        str
                RETURNS CITATION.

        """
    return (
        r" \line \line 1. Nanni, A., Titus-McQuillan, J., Moskalenko, O., "
        r"Pardo-Palacios, F., Liu, Z., Conesa, A., Rogers, R. L., & McIntyre, "
        "L. M. (2021). The evolution of splicing: transcriptome complexity "
        "and transcript distances implemented in TranD. \i bioRxiv\i0, "
        "2021.2009.2028.462251. https://doi.org/10.1101/2021.09.28.462251. "
        r"https://github.com/McIntyre-Lab/TranD."
    )


def endRtf(outFile):
    """
        Support function used to output legend text as RTF file.

        Parameters
        ----------
        outFile : FILE
                OUTPUT RTF FILE.

        Returns
        -------
        None.

        """
    outFile.write(
        r" Transcriptome analyses performed by TranD [1]. {}"
        r" \line \line \line \i Disclaimer:  While automated captions of "
        r"TranD have been carefully constructed, users are advised to verify "
        r"caption contents before use and report any errors to the TranD "
        r"github.\i0 ".format(
            getCitation()
        )
    )
    outFile.write(r"}\n\x00")     

#run the program
def main():
    ignoreDictionary = buildIgnoreDct(args)
    
    # check that there are >1 columns selected
    num_ignored = 0
    for ignore_flag in ignoreDictionary.values():
            if ignore_flag[0]:
                    num_ignored += 1


    if num_ignored > 4:
            exit("ERROR: There must be 2+ types of AS to plot. You can ignore"
                 " up to 4 types of AS.")
            
    # input csv to dataframe        
    inputDf = pd.read_csv(args.indir)

    trand.io.prepare_outdir(args.outdir, args.force)

    # retrieve name of input file
    input_file_name = os.path.splitext(os.path.basename(args.indir))[0]

    # prefix and ouput file prep
    output_file_name = "custom_{}{}".format(buildPrefix(args, ignoreDictionary), input_file_name)
    output_rtf_file = "{}/{}.rtf".format(args.outdir, output_file_name)
    
    #plot and output
    plotCustomPlot(inputDf, output_rtf_file, ignoreDictionary)
    plt.savefig("{}/{}.png".format(args.outdir, output_file_name), dpi=600, format="png")
    
    return 'KSB'


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

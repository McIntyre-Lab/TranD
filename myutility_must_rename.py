# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Thu Jan 26 12:43:56 2023

@author: k.bankole
"""

import trand.plot_functions as PF
#import plot_functions as PF
import trand.io
import argparse
import numpy as np
import pandas as pd
from IPython.display import display
from sys import exit

# think about 1 vs 2 GTF
# set arguments


def getOptions():
    # Parse command line arguments
    # Top part is adelena's description, change if neccessary
    parser = argparse.ArgumentParser(description="Make Upset plot for number of transcript pairs with each kind of AS"
                                     "for a single transcriptome with added box plots"
                                     "of the number of nt different and the proportion of nt different"
                                     "Make TranD plots from a GTF? file containing analysis of AS (--indir)"
                                     "Inludes the option to remove 3/5' variations (--remove-3-5)"
                                     "and/or No Shared NT from the plot (--rem-nosnt)")

    # Input data
    parser.add_argument(
        "-i",
        "-indir",
        dest="indir",
        required=True,
        help="Location of csv file containing transcriptome AS analysis")

    parser.add_argument(
        "-r35",
        "-remove-3-5",
        dest="rem35",
        action="store_true",
        help="Removes the 3' and 5' AS variations from the plot")

    parser.add_argument(
        "-rNSNT",
        "-remove-NoSharedNT",
        dest="remNSNT",
        action="store_true",
        help="Removes the transcript comparisons with no shared nucleotides")

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
        "-x",
        "--prefix",
        dest="outPrefix",
        required=True,
        help=(
            "Output prefix for plots. "
            # change this haha
            "Default: 'pre'"
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


# move to plot_functions?
# need to edit to add the choice to remove 53, noshared nt, or both, or neither!
# taking function from dadada, note modifications
# new name
def plot_pair_AS_upset_nt_box(td_data, legendOut, drop_5_3=False, drop_no_shared_nt=False):
    """
    Upset plot for the number of transcript pairs with each kind of AS
            for a single transcriptome with added box plots
            of the number of nt different and the proportion of nt different
    """
    pairAS = td_data[
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
    # Set No Shared NT pairs to have np.nan nucleotide difference
    pairAS["num_nt_diff"] = np.where(
        pairAS["flag_no_shared_nt"] == 1,
        np.nan,
        pairAS["num_nt_diff"]
    )
    # Make AS flags boolean values and rename
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

    if drop_5_3:
        # columns
        AScolsSubset = [
            "Alt. Donor/Acceptor",
            "Alt. Exon",
            "Intron Retention",
            "No Shared NT",
        ]
        # plot
        PF.plot_upset(
            pairAS.set_index(AScolsSubset),
            "",
            [
                "# NT Different",
                "Proportion\nNT Different",
            ],
        )
        # legend
        legendText = (
            "Number of transcript pairs with the specified types of alternative "
            "splicing indicated by the black dots below the histogram of transcript "
            "pair counts (n = {} transcript pairs, in {} multi-transcript genes). "
            "Box plots represent the number (blue) and proportion (orange) of "
            "nucleotides different between the pair. Transcript pairs with 5'/3' "
            "transcript length variation and/or nonoverlapping coordinates "
            "(and no IR or altnerative exon/donor/acceptor) are "
            "counted in the first column with no black dots.".format(
                len(pairAS),
                pairAS["gene_id"].nunique()
            )
        )
    elif drop_no_shared_nt:
        # columns
        AScolsSubset = [
            "5' Variation"
            "3' Variation"
            "Alt. Donor/Acceptor",
            "Alt. Exon",
            "Intron Retention",
            "No Shared NT",
        ]
        # plot
    elif drop_5_3 & drop_no_shared_nt:
        # columns
        AScolsSubset = [
            "Alt. Donor/Acceptor",
            "Alt. Exon",
            "Intron Retention",
        ]
        # plot
    else:
        # columns
        AScols = [
            "5' Variation",
            "3' Variation",
            "Alt. Donor/Acceptor",
            "Alt. Exon",
            "Intron Retention",
            "No Shared NT",
        ]
        # plot
        PF.plot_upset(
            pairAS.set_index(AScols),
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
            "nucleotides different between the pair. Transcript pairs with "
            "\"No Shared NT\" have nonoverlapping coordinates.".format(
                len(pairAS),
                pairAS["gene_id"].nunique()
            )
        )
    with open(legendOut, "w") as outFile:
        PF.start_rtf(outFile)
        outFile.write(
            r"\b Figure. Alternative splicing between pairs of transcripts \b0"
            r" \line {}".format(
                legendText
            )
        )
        PF.end_rtf(outFile)

        return pairAS


def check_args(args):
    del35 = False
    delNSNT = False

    del35 = args.rem35
    delNSNT = args.remNSNT

    return del35, delNSNT


# dont forget to validate inputs

def main():

    del35, delNSNT = check_args(args)

    trand.io.prepare_outdir(args.outdir, args.force)

    inputDf = pd.read_csv(args.indir)

    if del35:
        # remove 3' and 5' variations
        print("lolz: " + str(del35))
    elif delNSNT:
        # remove no shared nucleotide
        print("lmaoz: " + str(delNSNT))
    else:
        cool = "{} + agkonaognasgpo".format(args.outdir)
        print(cool)
        tempDf = plot_pair_AS_upset_nt_box(inputDf, args.outdir)

        print(tempDf)
    return 'KSB'


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:32:42 2023

@author: adalena/alison
"""

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Add variable to flag if small nt difference."
        )
    )

    # Input data
    parser.add_argument(
        "-d",
        "--distance-file",
        dest="inDist",
        required=True,
        help=(
            "Input file for TranD pairwise distance between the two sets of annotations."
        )
    )
    parser.add_argument(
        "-s",
        "--small-nt-difference",
        dest="smallDiff",
        required=False,
        type=int,
        help="Threshold for the number of nucleotides to use to define a 'small' difference in ERS transcript pairs."
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output pair classification file."
    )

    args = parser.parse_args()
    return args

def get_internal_nt_noOvlp(td_df):
    # Count number of internal (altnerative donor/acceptor or IR) nt differences in ERS
    #   Get the number of nt unique to T1/T2 in shared exon regions
    #   Subtract from this number the number of nt difference from a 5'/3' end length difference
    #   Only counting internal nt difference

    # Get difference of TSS
    #   Conditions:
    #       1. ERS, positive strand, 5' end difference, T1 has the longer 5' end
    #       2. ERS, positive strand, 5' end difference, T2 has the longer 5' end
    #       3. ERS, negative strand, 5' end difference, T1 has the longer 5' end
    #       4. ERS, negative strand, 5' end difference, T2 has the longer 5' end
    tssConditions = [
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_5_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[1] == td_df["EF_only_T1"].str.split("|").str[0].str.split(":").str[2]),
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_5_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[1] == td_df["EF_only_T2"].str.split("|").str[0].str.split(":").str[2]),
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_5_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="-")
           & (td_df["EF_ovlp"].str.split("|").str[-1].str.split(":").str[2] == td_df["EF_only_T1"].str.split("|").str[-1].str.split(":").str[1]),
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_5_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="-")
            & (td_df["EF_ovlp"].str.split("|").str[-1].str.split(":").str[2] == td_df["EF_only_T2"].str.split("|").str[-1].str.split(":").str[1]),

    ]
    tssChoices = [
        td_df["EF_only_T1"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["EF_only_T1"].str.split("|").str[0].str.split(":").str[1].astype(float),
        td_df["EF_only_T2"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["EF_only_T2"].str.split("|").str[0].str.split(":").str[1].astype(float),
        td_df["EF_only_T1"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["EF_only_T1"].str.split("|").str[-1].str.split(":").str[1].astype(float),
        td_df["EF_only_T2"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["EF_only_T2"].str.split("|").str[-1].str.split(":").str[1].astype(float),
    ]
    td_df["num_ERS_nt_noOvlp_TSS"] = np.select(tssConditions, tssChoices, np.nan)

    # Get difference of TTS
    #   Conditions:
    #       1. ERS, positive strand, 3' end difference, T1 has the longer 3' end
    #       2. ERS, positive strand, 3' end difference, T2 has the longer 3' end
    #       3. ERS, negative strand, 3' end difference, T1 has the longer 3' end
    #       4. ERS, negative strand, 3' end difference, T2 has the longer 3' end
    ttsConditions = [
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_3_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["EF_ovlp"].str.split("|").str[-1].str.split(":").str[2] == td_df["EF_only_T1"].str.split("|").str[-1].str.split(":").str[1]),
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_3_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["EF_ovlp"].str.split("|").str[-1].str.split(":").str[2] == td_df["EF_only_T2"].str.split("|").str[-1].str.split(":").str[1]),
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_3_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="-")
           & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[1] == td_df["EF_only_T1"].str.split("|").str[0].str.split(":").str[2]),
        (td_df["prop_ER_noOvlp"] == 0)
            & (td_df["flag_3_var"]==1)
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[3]=="-")
            & (td_df["EF_ovlp"].str.split("|").str[0].str.split(":").str[1] == td_df["EF_only_T2"].str.split("|").str[0].str.split(":").str[2]),

    ]
    ttsChoices = [
        td_df["EF_only_T1"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["EF_only_T1"].str.split("|").str[-1].str.split(":").str[1].astype(float),
        td_df["EF_only_T2"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["EF_only_T2"].str.split("|").str[-1].str.split(":").str[1].astype(float),
        td_df["EF_only_T1"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["EF_only_T1"].str.split("|").str[0].str.split(":").str[1].astype(float),
        td_df["EF_only_T2"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["EF_only_T2"].str.split("|").str[0].str.split(":").str[1].astype(float),
    ]
    td_df["num_ERS_nt_noOvlp_TTS"] = np.select(ttsConditions, ttsChoices, np.nan)

    # Get donor/acceptor/IR length difference
    td_df["num_ERS_nt_noOvlp_internal"] = np.where(
        td_df["prop_ER_noOvlp"] == 0,
        td_df["num_nt_noOvlp"].astype(int) - td_df["num_ERS_nt_noOvlp_TSS"].fillna(0) - td_df["num_ERS_nt_noOvlp_TTS"].fillna(0),
        np.nan
    )
    return td_df

def main():
    # TranD 2GTF distance of consolidated dataset1 vs. dataset2
    distDf = pd.read_csv(args.inDist, low_memory=False)

    # Get internal nt difference
    distDf2 = get_internal_nt_noOvlp(distDf)
    del(distDf)

    # Classify transcript pairs
    # If requested add classification that includes "small" ERS differences (provide a max # of nt in ERS_noIR)
    #rmpDf = distDf2[distDf2["flag_RMP"]==1].copy()
    #print(rmpDf)

    if args.smallDiff is not None:
        distDf2["flag_ERS_noIR"] = np.where(
            (distDf2["prop_ER_noOvlp"]==0) & (distDf2["flag_IR"]==0),
            1,
            0
        )
        distDf2["flag_ERS_wIR"] = np.where(
            (distDf2["prop_ER_noOvlp"]==0) & (distDf2["flag_IR"]==1),
            1,
            0
        )
        #distDf2["flag_RMP"] = 1
        #print(distDf2.columns.values)

        # Add variable that is FSM, ERS_noIR_small (ERS_noIR with < N nt internal difference), ERS_noIR_large (ERS_noIR with >= N nt internal difference), ERS_wIR, ERN (recip min that is not FSM/ERS), NRM (no reciprocal minimum match)
        compConditions = [
            distDf2["flag_FSM"] == 1,
            (distDf2["flag_ERS_noIR"] == 1) & (distDf2["num_ERS_nt_noOvlp_internal"] < args.smallDiff),
            (distDf2["flag_ERS_noIR"] == 1) & (distDf2["num_ERS_nt_noOvlp_internal"] >= args.smallDiff),
            distDf2["flag_ERS_wIR"]==1,
            distDf2["flag_RMP"] == 1,
            distDf2["flag_RMP"] == 0,
        ]
        compChoices = [
            "FSM",
            "ERS_noIR_small",
            "ERS_noIR_large",
            "ERS_wIR",
            "ERN",
            "NRM",
        ]        
    else:
        # Add variable that is FSM, ERS_noIR, ERS_wIR, ERN (recip min that is not FSM/ERS), NRM (no reciprocal minimum match)
        compConditions = [
            distDf2["flag_FSM"] == 1,
            distDf2["flag_ERS_noIR"] == 1,
            distDf2["flag_ERS_wIR"]==1,
            distDf2["flag_RMP"] == 1,
            distDf2["flag_RMP"] == 0,
        ]
        compChoices = [
            "FSM",
            "ERS_noIR",
            "ERS_wIR",
            "ERN",
            "NRM",
        ]
    distDf2["pair_classification"] = np.select(compConditions, compChoices, "oops")

    outDf = distDf2[["gene_id", "transcript_1", "transcript_2", "flag_FSM", "flag_RMP", "flag_ERS_noIR", "pair_classification"]]

    # Output transcript map file
    outDf.to_csv(args.outFile, index=False)
#    unionRMP.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/hg38_RefSeq_Ensembl_transcript_map.csv", index=False)

    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

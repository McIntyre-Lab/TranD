#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:32:42 2023

@author: adalena
"""

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Create a Transcript Map File from TranD output files."
        )
    )

    # Input data
    parser.add_argument(
        "-k",
        "--key-file",
        dest="inKey",
        required=True,
        help=(
            "Input file for union key file between the two sets of annotations."
        )
    )
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
        "-e",
        "--erg-2gtf",
        dest="erGroup",
        required=False,
        help=(
            "Input file of exon region group (ERG) information based on the pairwise distance of the two sets of annotations."
        )
    )
    parser.add_argument(
        "-1",
        "--name1",
        dest="inName1",
        required=True,
        help=(
            "Input annotation 1 name."
        )
    )
    parser.add_argument(
        "--prefix1",
        dest="inPrefix1",
        required=False,
        help="Prefix added to consolidated transcripts for dataset 1."
    )
    parser.add_argument(
        "-2",
        "--name2",
        dest="inName2",
        required=True,
        help=(
            "Input annotation 2 name."
        )
    )
    parser.add_argument(
        "--prefix2",
        dest="inPrefix2",
        required=False,
        help="Prefix added to consolidated transcripts for dataset 2."
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
        help="Output Transcript Map file."
    )

    args = parser.parse_args()
    return args

def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
    # Split variable by some character like '|' or ',' and keep all other values the same
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def get_internal_nt_diff(td_df):
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
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_5_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[1] == td_df["fragment_T1_only"].str.split("|").str[0].str.split(":").str[2]),
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_5_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[1] == td_df["fragment_T2_only"].str.split("|").str[0].str.split(":").str[2]),
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_5_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="-")
           & (td_df["fragment_shared"].str.split("|").str[-1].str.split(":").str[2] == td_df["fragment_T1_only"].str.split("|").str[-1].str.split(":").str[1]),
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_5_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="-")
            & (td_df["fragment_shared"].str.split("|").str[-1].str.split(":").str[2] == td_df["fragment_T2_only"].str.split("|").str[-1].str.split(":").str[1]),

    ]
    tssChoices = [
        td_df["fragment_T1_only"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["fragment_T1_only"].str.split("|").str[0].str.split(":").str[1].astype(float),
        td_df["fragment_T2_only"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["fragment_T2_only"].str.split("|").str[0].str.split(":").str[1].astype(float),
        td_df["fragment_T1_only"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["fragment_T1_only"].str.split("|").str[-1].str.split(":").str[1].astype(float),
        td_df["fragment_T2_only"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["fragment_T2_only"].str.split("|").str[-1].str.split(":").str[1].astype(float),
    ]
    td_df["num_ERS_nt_diff_TSS"] = np.select(tssConditions, tssChoices, np.nan)

    # Get difference of TTS
    #   Conditions:
    #       1. ERS, positive strand, 3' end difference, T1 has the longer 3' end
    #       2. ERS, positive strand, 3' end difference, T2 has the longer 3' end
    #       3. ERS, negative strand, 3' end difference, T1 has the longer 3' end
    #       4. ERS, negative strand, 3' end difference, T2 has the longer 3' end
    ttsConditions = [
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_3_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["fragment_shared"].str.split("|").str[-1].str.split(":").str[2] == td_df["fragment_T1_only"].str.split("|").str[-1].str.split(":").str[1]),
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_3_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="+")
            & (td_df["fragment_shared"].str.split("|").str[-1].str.split(":").str[2] == td_df["fragment_T2_only"].str.split("|").str[-1].str.split(":").str[1]),
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_3_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="-")
           & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[1] == td_df["fragment_T1_only"].str.split("|").str[0].str.split(":").str[2]),
        (td_df["prop_ER_diff"] == 0)
            & (td_df["flag_3_variation"]==1)
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[3]=="-")
            & (td_df["fragment_shared"].str.split("|").str[0].str.split(":").str[1] == td_df["fragment_T2_only"].str.split("|").str[0].str.split(":").str[2]),

    ]
    ttsChoices = [
        td_df["fragment_T1_only"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["fragment_T1_only"].str.split("|").str[-1].str.split(":").str[1].astype(float),
        td_df["fragment_T2_only"].str.split("|").str[-1].str.split(":").str[2].astype(float) - td_df["fragment_T2_only"].str.split("|").str[-1].str.split(":").str[1].astype(float),
        td_df["fragment_T1_only"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["fragment_T1_only"].str.split("|").str[0].str.split(":").str[1].astype(float),
        td_df["fragment_T2_only"].str.split("|").str[0].str.split(":").str[2].astype(float) - td_df["fragment_T2_only"].str.split("|").str[0].str.split(":").str[1].astype(float),
    ]
    td_df["num_ERS_nt_diff_TTS"] = np.select(ttsConditions, ttsChoices, np.nan)

    # Get donor/acceptor/IR length difference
    td_df["num_ERS_nt_diff_internal"] = np.where(
        td_df["prop_ER_diff"] == 0,
        td_df["num_nt_diff"].astype(int) - td_df["num_ERS_nt_diff_TSS"].fillna(0) - td_df["num_ERS_nt_diff_TTS"].fillna(0),
        np.nan
    )
    return td_df

def main():
    # Get annotation names
    name1 = args.inName1
    name2 = args.inName2
#    name1 = "RefSeq"
#    name2 = "Ensembl"

    # Get prefixes used for consolidation of each annotation
    prefix1 = args.inPrefix1
    prefix2 = args.inPrefix2
#    prefix1 = "hg38_refseq_2_ensembl"
#    prefix2 = "hg38_ensembl"

    # Union key file
    unionDf = pd.read_csv(args.inKey)
#    unionDf = pd.read_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/hg38_refseq_ensembl_union_ref_keyFile.csv")
    
    # Flag UJC in one or both annotations
    unionDf["flag_in_"+name1+"_only"] = np.where(
        (~unionDf[name1+"_transcript_id"].isna()) & (unionDf[name2+"_transcript_id"].isna()),
        1,
        0
    )
    unionDf["flag_in_"+name2+"_only"] = np.where(
        (~unionDf[name2+"_transcript_id"].isna()) & (unionDf[name1+"_transcript_id"].isna()),
        1,
        0
    )
    unionDf["flag_in_both"] = np.where(
        (~unionDf[name1+"_transcript_id"].isna()) & (~unionDf[name2+"_transcript_id"].isna()),
        1,
        0
    )
    if len(unionDf.fillna("")[unionDf.fillna("")[name2+"_uniq_jxn_id"].str.contains("\|")]) > 0:
        print("WARNING: There is a union UJC that covers more than one {} UJC -\n{}".format(name2,unionDf.fillna("")[unionDf.fillna("")[name2+"_uniq_jxn_id"].str.contains("\|")].to_string(index=False)))
    if len(unionDf.fillna("")[unionDf.fillna("")[name1+"_uniq_jxn_id"].str.contains("\|")]) > 0:
        print("WARNING: There is a union UJC that covers more than one {} UJC -\n{}".format(name1,unionDf.fillna("")[unionDf.fillna("")[name1+"_uniq_jxn_id"].str.contains("\|")].to_string(index=False)))

    # Flag genes that are only in one annotation or the other
    geneUnion = unionDf.groupby("gene_id")[["flag_in_"+name1+"_only", "flag_in_"+name2+"_only", "flag_in_both"]].max()
    geneChoices = [name1, name2, "both"]
    geneConditions = [
        (geneUnion.sum(axis=1)==1) & (geneUnion["flag_in_"+name1+"_only"]==1),
        (geneUnion.sum(axis=1)==1) & (geneUnion["flag_in_"+name2+"_only"]==1),
        ((geneUnion.sum(axis=1)==1) & (geneUnion["flag_in_both"]==1)) | (geneUnion.sum(axis=1)>1),
    ]
    geneUnion["gene_annotation"] = np.select(geneConditions, geneChoices, "oops")
    unionMergeDf = pd.merge(
        unionDf,
        geneUnion.reset_index()[["gene_id", "gene_annotation"]],
        how="outer",
        on="gene_id",
    )

    # TranD 2GTF distance of consolidated dataset1 vs. dataset2
    distDf = pd.read_csv(args.inDist, low_memory=False)
#    distDf = pd.read_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/TranD_consol_refseq_vs_ensembl/all_chrom_pairwise_transcript_distance.csv", low_memory=False)

    # Get internal nt difference
    distDf2 = get_internal_nt_diff(distDf)
    del(distDf)

    # Drop all distance rows that are not minimum pairs and remove suffix from TranD 2 GTF pairwise
    distDf3 = distDf2[distDf2["flag_min_match_"+prefix1]+distDf2["flag_min_match_"+prefix2]>0].copy()
    del(distDf2)
    distDf3[name1+"_uniq_jxn_id"] = distDf3["transcript_1"].str[:-len(prefix1)-1]
    distDf3[name2+"_uniq_jxn_id"] = distDf3["transcript_2"].str[:-len(prefix2)-1]        

    # Get all pairs that are reciprocal minimum and extras
    rmpDf = distDf3[distDf3["flag_recip_min_match"]==1].copy()
    rmpDf["flag_ERS_noIR"] = np.where(
        (rmpDf["prop_ER_diff"]==0) & (rmpDf["flag_IR"]==0),
        1,
        0
    )
    rmpDf["flag_ERS_withIR"] = np.where(
        (rmpDf["prop_ER_diff"]==0) & (rmpDf["flag_IR"]==1),
        1,
        0
    )
    rmpDf["flag_RMP"] = 1

    # !!! TODO: merge in minimum distance of "extras" that do not have RMP
#    extraDf = distDf3[distDf3["flag_recip_min_match"]!=1].copy()
    
    # Merge minimum distance pairs with union UJC in both annotations
    unionRMP = pd.merge(
        rmpDf,
        unionMergeDf[unionMergeDf["flag_in_both"]==1].drop(columns=[name2+"_transcript_id", "flag_in_"+name2+"_only", "flag_in_"+name1+"_only", "flag_in_both"]),
        how="outer",
        on=["gene_id", name1+"_uniq_jxn_id", name2+"_uniq_jxn_id"],
        indicator="merge_check",
        validate="1:1"
    )
#    unionRMP["merge_check"].value_counts()
    # For hg38 refseq vs. ensembl example
    # left_only      56505
    # both           10290
    # right_only        15
    if unionRMP["merge_check"].value_counts()["right_only"] > 0:
        print("WARNING: There are pairs of UJC_id that are not found in the RMP file:\n{}".format(unionRMP[unionRMP["merge_check"]=="right_only"].to_csv(index=False)))
    
    pairUnionRMP = unionRMP[unionRMP["merge_check"]=="both"]
    noPairUnionRMP = unionRMP[unionRMP["merge_check"]=="left_only"][[c for c in rmpDf.columns if c not in ["transcript_1", "transcript_2"]]]
    
    # Merge RMP without a union UJC pair with the RMP that are only in species1 and then those that are only in species2
    unionRMP1 = pd.merge(
        unionMergeDf[unionMergeDf["flag_in_"+name1+"_only"]==1].drop(columns=[name2+"_uniq_jxn_id", name2+"_transcript_id", "flag_in_"+name2+"_only", "flag_in_"+name1+"_only", "flag_in_both"]),
        noPairUnionRMP,
        how="outer",
        on=["gene_id", name1+"_uniq_jxn_id"],
        indicator="merge_check",
        validate="1:1"
    )
    if unionRMP1["merge_check"].value_counts()["right_only"] > 0:
        print("WARNING: There are UJC_id for {} in the RMP file that are not found in the UJC file:\n{}".format(name1, unionRMP1[unionRMP1["merge_check"]=="right_only"].to_csv(index=False)))
    
    unionRMP1Present = unionRMP1[unionRMP1["merge_check"]=="both"].drop(columns=["merge_check"])
    unionRMP1No = unionRMP1[unionRMP1["merge_check"]=="left_only"].drop(columns=["merge_check"])
    unionRMP2 = pd.merge(
        unionMergeDf[unionMergeDf["flag_in_"+name2+"_only"]==1].drop(columns=[name1+"_uniq_jxn_id", name1+"_transcript_id", "flag_in_"+name2+"_only", "flag_in_"+name1+"_only", "flag_in_both"]),
        unionRMP1Present.rename(columns={"union_uniq_jxn_id": name1+"_union_uniq_jxn_id"}),
        how="outer",
        on=["gene_id", "gene_annotation", name2+"_uniq_jxn_id"],
        indicator="merge_check",
        validate="1:1"
    )
    if unionRMP2["merge_check"].value_counts()["right_only"] > 0:
        print("WARNING: There are UJC_id for {} in the RMP file that are not found in the UJC file:\n{}".format(name2, unionRMP2[unionRMP2["merge_check"]=="right_only"].to_csv(index=False)))
    
    # Combine merged RMP union UJC pairs, merged individual annotation RMP UJC pairs, and inidvidual annotation unmerged UJC
    pairUnionRMP = pairUnionRMP.rename(columns={"union_uniq_jxn_id": name1+"_union_uniq_jxn_id"})
    pairUnionRMP[name2+"_union_uniq_jxn_id"] = pairUnionRMP[name1+"_union_uniq_jxn_id"]
    unionRMP1No = unionRMP1No.rename(columns={"union_uniq_jxn_id": name1+"_union_uniq_jxn_id"})
    unionRMP2 = unionRMP2.rename(columns={"union_uniq_jxn_id": name2+"_union_uniq_jxn_id"})
    unionRMP = pd.concat([pairUnionRMP, unionRMP1No, unionRMP2], ignore_index=True).drop(columns=["merge_check"])
    
    # Fill in na values where necessary
    unionRMP[[c for c in unionRMP.columns if "flag_" in c]] = unionRMP[[c for c in unionRMP.columns if "flag_" in c]].fillna(0)

    # Classify transcript pairs
    # If requested add classification that includes "small" ERS differences (provide a max # of nt in ERS_noIR)
    if args.smallDiff is not None:
        # Add variable that is FSM, ERS_noIR_small (ERS_noIR with < N nt internal difference), ERS_noIR_large (ERS_noIR with >= N nt internal difference), ERS_wIR, ERN (recip min that is not FSM/ERS), NRM (no reciprocal minimum match)
        compConditions = [
            unionRMP["flag_FSM"] == 1,
            (unionRMP["flag_ERS_noIR"] == 1) & (unionRMP["num_ERS_nt_diff_internal"] < args.smallDiff),
            (unionRMP["flag_ERS_noIR"] == 1) & (unionRMP["num_ERS_nt_diff_internal"] >= args.smallDiff),
            unionRMP["flag_ERS_withIR"]==1,
            unionRMP["flag_RMP"] == 1,
            (unionRMP["flag_RMP"] == 0) & (unionRMP[name1+"_transcript_id"].isna()),
            (unionRMP["flag_RMP"] == 0) & (unionRMP[name2+"_transcript_id"].isna()),    
        ]
        compChoices = [
            "FSM",
            "ERS_noIR_small",
            "ERS_noIR_large",
            "ERS_wIR",
            "ERN",
            "NRM_"+name2,
            "NRM_"+name1,
        ]        
    else:
        # Add variable that is FSM, ERS_noIR, ERS_wIR, ERN (recip min that is not FSM/ERS), NRM (no reciprocal minimum match)
        compConditions = [
            unionRMP["flag_FSM"] == 1,
            unionRMP["flag_ERS_noIR"] == 1,
            unionRMP["flag_ERS_withIR"]==1,
            unionRMP["flag_RMP"] == 1,
            (unionRMP["flag_RMP"] == 0) & (unionRMP[name1+"_transcript_id"].isna()),
            (unionRMP["flag_RMP"] == 0) & (unionRMP[name2+"_transcript_id"].isna()),    
        ]
        compChoices = [
            "FSM",
            "ERS_noIR",
            "ERS_wIR",
            "ERN",
            "NRM_"+name2,
            "NRM_"+name1,
        ]
    unionRMP["comparison_type"] = np.select(compConditions, compChoices, "oops")
    
    # Fill in transcript_1 and transcript_2 values
    unionRMP["transcript_1"] = np.where(
        unionRMP["transcript_1"].isna(),
        unionRMP[name1+"_uniq_jxn_id"]+"_"+prefix1,
        unionRMP["transcript_1"]
    )
    unionRMP["transcript_2"] = np.where(
        unionRMP["transcript_2"].isna(),
        unionRMP[name2+"_uniq_jxn_id"]+"_"+prefix2,
        unionRMP["transcript_2"]
    )

    # Get ER group output if provided and merge with distance
    if args.erGroup is not None:
        ergDf = split_column_by_sep(
            pd.read_csv(args.erGroup), col_name="xscripts", sep="|")
        mergeRmp1 = pd.merge(
            ergDf[["xscripts", "contains_which_gtf"]].rename(columns={"contains_which_gtf": "contains_which_gtf_"+name1}),
            unionRMP,
            how="outer",
            left_on="xscripts",
            right_on="transcript_1",
            validate="1:m",
            indicator="merge_check"
        )
        mergeRmp2 = pd.merge(
            mergeRmp1[mergeRmp1["merge_check"]=="left_only"][["xscripts", "contains_which_gtf_"+name1]].rename(columns={"contains_which_gtf_"+name1: "contains_which_gtf_"+name2}),
            mergeRmp1[mergeRmp1["merge_check"]!="left_only"][list(unionRMP.columns) + ["contains_which_gtf_"+name1]],
            how="outer",
            left_on="xscripts",
            right_on="transcript_2",
            validate="1:m",
            indicator="merge_check"
        )
        unionRMP = mergeRmp2[mergeRmp2["merge_check"]!="left_only"].drop(columns=["merge_check", "xscripts"])
        ergChoices = [name1 , name2 , "both"]
        ergConditions1 = [
            unionRMP["contains_which_gtf_"+name1]==1,
            unionRMP["contains_which_gtf_"+name1]==2,
            unionRMP["contains_which_gtf_"+name1]==3,            
        ]
        ergConditions2 = [
            unionRMP["contains_which_gtf_"+name2]==1,
            unionRMP["contains_which_gtf_"+name2]==2,
            unionRMP["contains_which_gtf_"+name2]==3,            
        ]
        unionRMP["erg_content_"+name1] = np.select(ergConditions1, ergChoices, "oops")
        del(unionRMP["contains_which_gtf_"+name1])
        unionRMP["erg_content_"+name2] = np.select(ergConditions2, ergChoices, "oops")
        del(unionRMP["contains_which_gtf_"+name2])

        # Fill in ERG for transcripts in genes only in one of the two annotations
        unionRMP["erg_content_"+name1] = np.select(
            [unionRMP["erg_content_"+name1] != "oops",
             (unionRMP["erg_content_"+name1]=="oops") & (unionRMP["transcript_2"].isna()),
             (unionRMP["erg_content_"+name1]=="oops") & (unionRMP["transcript_1"].isna())],
            [unionRMP["erg_content_"+name1], name1, name2],
            "oops"
        )
        unionRMP["erg_content_"+name2] = np.select(
            [unionRMP["erg_content_"+name2] != "oops",
             (unionRMP["erg_content_"+name2]=="oops") & (unionRMP["transcript_2"].isna()),
             (unionRMP["erg_content_"+name2]=="oops") & (unionRMP["transcript_1"].isna())],
            [unionRMP["erg_content_"+name2], name1, name2],
            "oops"
        )        

    # Output transcript map file
    unionRMP.to_csv(args.outFile, index=False)
#    unionRMP.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/human_analysis/hg38_RefSeq_Ensembl_transcript_map.csv", index=False)

    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
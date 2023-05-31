#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 12:49:29 2021

@author: adalena.nanni
"""

import numpy as np
from loguru import logger


def get_md_cols(name1, name2):
    md_df_cols = [
            'gene_id',
            'transcript_1',
            'transcript_2',
            'num_jxn_only_T1',
            'num_jxn_only_T2',
            'num_jxn_ovlp',
            'prop_jxn_noOvlp',
            'prop_jxn_ovlp',
            'jxn_string_T1',
            'jxn_string_T2',
            'jxn_only_T1',
            'jxn_only_T2',
            'jxn_same',
            'num_ER_only_T1',
            'num_ER_only_T2',
            'num_ER_ovlp',
            'prop_ER_noOvlp',
            'prop_ER_ovlp',
            'ER_only_T1',
            'ER_only_T2',
            'ER_ovlp',
            'num_EF_only_T1',
            'num_EF_only_T2',
            'num_EF_ovlp',
            'prop_EF_noOvlp',
            'prop_EF_ovlp',
            'EF_only_T1',
            'EF_only_T2',
            'EF_ovlp',
            'num_exon_only_T1',
            'num_exon_only_T2',
            'num_exon_ovlp',
            'num_IR_EF_T1',
            'num_IR_EF_T2',
            'IR_EF_T1',
            'IR_EF_T2',
            'num_nt_ovlp',
            'num_nt_only_T1',
            'num_nt_only_T2',
            'num_nt_noOvlp',
            'total_nt',
            'prop_nt_noOvlp',
            'prop_nt_ovlp',
            'num_nt_only_T1_in_ovlpER',
            'num_nt_only_T2_in_ovlpER',
            'num_nt_ovlp_in_ovlpER',
            'total_nt_in_ovlpER',
            'prop_nt_noOvlp_in_ovlpER',
            'prop_nt_ovlp_in_ovlpER',
            'num_nt_only_T1_in_uniqER',
            'num_nt_only_T2_in_uniqER',
            'flag_FSM',
            'flag_T1_ISM_of_T2',
            'flag_T2_ISM_of_T1',
            'flag_IR',
            'flag_5_var',
            'flag_3_var',
            'flag_alt_DA',
            'flag_alt_exon',
            'flag_no_ovlp_nt',
            'num_transcript_inGene_'+name1,
            'num_transcript_inGene_'+name2,
            'flag_' + name1 + '_greater',
            'flag_' + name2 + '_greater',
            'flag_match',
            'transcript_inGene',
            'min_match_' + name1,
            'flag_min_match_' + name1,
            'min_match_' + name2,
            'flag_min_match_' + name2,
            'flag_RMP',
            'num_RMP_inGene',
            'flag_' + name1 + '_tie',
            'flag_' + name2 + '_tie',
            name1 + '_distance_tie',
            name2 + '_distance_tie',
            'flag_identical_RMP',
            'flag_FSM_RMP',
            'flag_ERS_RMP',
            'flag_ERS_noIR_RMP',
            'flag_ERS_wIR_RMP',
            'num_identical_RMP_inGene',
            'num_FSM_RMP_inGene',
            'num_ERS_RMP_inGene',
            'num_ERS_noIR_RMP_inGene',
            'num_ERS_wIR_RMP_inGene',
            'prop_identical_RMP_inGene',
            'prop_FSM_RMP_inGene',
            'prop_ERS_RMP_inGene',
            'prop_ERS_noIR_RMP_inGene',
            'prop_ERS_wIR_RMP_inGene',
            'flag_no_ovlp_nt_RMP',
            'flag_alt_exon_RMP',
            'flag_alt_DA_RMP',
            'flag_IR_RMP',
            'flag_5_var_RMP',
            'flag_3_var_RMP',
            'recip_min_pair_inGene',
            'flag_min_match_' + name1 + '_w_tie',
            'flag_min_match_' + name2 + '_w_tie',
            'flag_RMP_w_tie',
            'num_RMP_inGene_w_tie',
            'flag_identical_RMP_w_tie',
            'flag_FSM_RMP_w_tie',
            'flag_ERS_RMP_w_tie',
            'flag_ERS_noIR_RMP_w_tie',
            'flag_ERS_wIR_RMP_w_tie',
            'num_identical_RMP_inGene_w_tie',
            'num_FSM_RMP_inGene_w_tie',
            'num_ERS_RMP_inGene_w_tie',
            'num_ERS_noIR_RMP_inGene_w_tie',
            'prop_identical_RMP_inGene_w_tie',
            'prop_FSM_RMP_inGene_w_tie',
            'prop_ERS_RMP_inGene_w_tie',
            'prop_ERS_noIR_RMP_inGene_w_tie',
            'num_ERS_wIR_RMP_inGene_w_tie',
            'prop_ERS_wIR_RMP_inGene_w_tie',
            'recip_min_pair_inGene_w_tie',
        ]
    return md_df_cols

def get_num_inGene_w_tie(df, flag_name):
    df2 = df.copy()
    df2["name_1_flag_w_tie"] = np.where(
           df2[flag_name] == 1,
           df2["transcript_1"],
           np.nan
    )
    df2["name_2_flag_w_tie"] = np.where(
           df2[flag_name] == 1,
           df2["transcript_2"],
           np.nan
    )
    df2["num_name_1_flag_w_tie"] = (
        df2.groupby("gene_id")[
            "name_1_flag_w_tie"
        ].transform("nunique")
    )
    df2["num_name_2_flag_w_tie"] = (
        df2.groupby("gene_id")[
            "name_2_flag_w_tie"
        ].transform("nunique")
    )
    df2["num_identical_RMP_inGene_w_tie"] = df2[[
            "num_name_1_flag_w_tie",
            "num_name_2_flag_w_tie"
        ]].min(axis=1)
    return df2["num_identical_RMP_inGene_w_tie"]

def identify_min_pair(td_data, out_pairs, name1, name2):

    # Get number of transcripts per gene
    td_data["num_transcript_inGene_" + name1] = td_data.groupby("gene_id")[
        "transcript_1"
    ].transform("nunique")
    td_data["num_transcript_inGene_" + name2] = td_data.groupby("gene_id")[
        "transcript_2"
    ].transform("nunique")

    # Flag groups based on transcripts per gene
    td_data["flag_" + name1 + "_greater"] = np.where(
        (
            td_data["num_transcript_inGene_" + name1]
            > td_data["num_transcript_inGene_" + name2]
        )
        & (td_data["num_transcript_inGene_" + name1] > 0)
        & (td_data["num_transcript_inGene_" + name2] > 0),
        1,
        0,
    )
    td_data["flag_" + name2 + "_greater"] = np.where(
        (
            td_data["num_transcript_inGene_" + name1]
            < td_data["num_transcript_inGene_" + name2]
        )
        & (td_data["num_transcript_inGene_" + name1] > 0)
        & (td_data["num_transcript_inGene_" + name2] > 0),
        1,
        0,
    )
    td_data["flag_match"] = np.where(
        (
            td_data["num_transcript_inGene_" + name1]
            == td_data["num_transcript_inGene_" + name2]
        ),
        1,
        0,
    )
    conditionsGene = [
        td_data["flag_" + name1 + "_greater"] == 1,
        td_data["flag_" + name2 + "_greater"] == 1,
        td_data["flag_match"] == 1,
    ]
    choicesGene = [name1 + "_greater", name2 + "_greater", "match"]
    td_data["transcript_inGene"] = np.select(conditionsGene, choicesGene, "missing")
    if len(td_data[td_data["transcript_inGene"] == "missing"]) > 0:
        logger.error(
            "Unexpected variable assignment for gene transcript count comparison"
        )

    # Flag the minimum pair for each transcript
    # First sort (ascending) by proportion of junctions different,
    #   proportion of ER different, and proportion of nt different
    # Also sort at the end by transcript_1 and transcrpit_2 IDs to ensure the
    #   same transcript is chosen in the case of ties
    td_data = td_data.sort_values(
        [
            "prop_ER_noOvlp",
            "prop_jxn_noOvlp",
            "prop_nt_noOvlp",
            "transcript_1",
            "transcript_2",
        ],
        ascending=[True, True, True, False, False],
    )
    td_data["min_match_" + name1] = td_data.groupby("transcript_1")[
        "transcript_2"
    ].transform("first")
    td_data["flag_min_match_" + name1] = np.where(
        td_data["min_match_" + name1] == td_data["transcript_2"], 1, 0
    )
    td_data["min_match_" + name2] = td_data.groupby("transcript_2")[
        "transcript_1"
    ].transform("first")
    td_data["flag_min_match_" + name2] = np.where(
        td_data["min_match_" + name2] == td_data["transcript_1"], 1, 0
    )
    td_data["flag_RMP"] = np.where(
        (td_data["flag_min_match_" + name1] == 1)
        & (td_data["flag_min_match_" + name2] == 1),
        1,
        0,
    )
    td_data["num_RMP_inGene"] = td_data.groupby("gene_id")[
        "flag_RMP"
    ].transform("sum")

    # Flag pairs that are ties in between at least 2 transcirpts in the other dataset
    td_data["flag_" + name1 + "_tie"] = (
        td_data[
            [
                "transcript_1",
                "prop_jxn_noOvlp",
                "prop_ER_noOvlp",
                "prop_nt_noOvlp",
                "num_nt_noOvlp",
            ]
        ]
        .duplicated(keep=False)
        .astype(int)
    )
    td_data["flag_" + name2 + "_tie"] = (
        td_data[
            [
                "transcript_2",
                "prop_jxn_noOvlp",
                "prop_ER_noOvlp",
                "prop_nt_noOvlp",
                "num_nt_noOvlp",
            ]
        ]
        .duplicated(keep=False)
        .astype(int)
    )

    # List the other transcripts that the distances are tied with
    td_data[name1 + "_distance_tie"] = td_data.groupby(
            [
                "transcript_1",
                "prop_jxn_noOvlp",
                "prop_ER_noOvlp",
                "prop_nt_noOvlp",
                "num_nt_noOvlp",
            ])["transcript_2"].transform(lambda x: "|".join(x))
    td_data[name2 + "_distance_tie"] = td_data.groupby(
            [
                "transcript_2",
                "prop_jxn_noOvlp",
                "prop_ER_noOvlp",
                "prop_nt_noOvlp",
                "num_nt_noOvlp",
            ])["transcript_1"].transform(lambda x: "|".join(x))
        
    # Flag FSM (all junctions matching) and ERS (all exonic regions shared) min matches
    td_data["flag_identical_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["prop_nt_noOvlp"] == 0),
        1,
        0
    )
    td_data["flag_FSM_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["prop_jxn_noOvlp"] == 0),
        1,
        0,
    )
    td_data["flag_ERS_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["prop_ER_noOvlp"] == 0),
        1,
        0
    )
    td_data["flag_ERS_noIR_RMP"] = np.where(
        (td_data["flag_RMP"] == 1)
        & (td_data["prop_ER_noOvlp"] == 0)
        & (td_data["num_IR_EF_T1"] + td_data["num_IR_EF_T2"] == 0),
        1,
        0,
    )
    td_data["flag_ERS_wIR_RMP"] = np.where(
        (td_data["flag_RMP"] == 1)
        & (td_data["prop_ER_noOvlp"] == 0)
        & (td_data["num_IR_EF_T1"] + td_data["num_IR_EF_T2"] > 0),
        1,
        0,
    )
    td_data["num_identical_RMP_inGene"] = td_data.groupby("gene_id")[
        "flag_identical_RMP"
    ].transform("sum")
    td_data["num_FSM_RMP_inGene"] = td_data.groupby("gene_id")[
        "flag_FSM_RMP"
    ].transform("sum")
    td_data["num_ERS_RMP_inGene"] = td_data.groupby("gene_id")[
        "flag_ERS_RMP"
    ].transform("sum")
    td_data["num_ERS_noIR_RMP_inGene"] = td_data.groupby("gene_id")[
        "flag_ERS_noIR_RMP"
    ].transform("sum")
    td_data["num_ERS_wIR_RMP_inGene"] = td_data.groupby("gene_id")[
        "flag_ERS_wIR_RMP"
    ].transform("sum")
    td_data["prop_identical_RMP_inGene"] = (
        td_data["num_identical_RMP_inGene"]
        / td_data["num_RMP_inGene"]
    )
    td_data["prop_FSM_RMP_inGene"] = (
        td_data["num_FSM_RMP_inGene"]
        / td_data["num_RMP_inGene"]
    )
    td_data["prop_ERS_RMP_inGene"] = (
        td_data["num_ERS_RMP_inGene"]
        / td_data["num_RMP_inGene"]
    )
    td_data["prop_ERS_noIR_RMP_inGene"] = (
        td_data["num_ERS_noIR_RMP_inGene"]
        / td_data["num_RMP_inGene"]
    )
    td_data["prop_ERS_wIR_RMP_inGene"] = (
        td_data["num_ERS_wIR_RMP_inGene"]
        / td_data["num_RMP_inGene"]
    )

    # Flag different types of splicing differences in min matches
    #   (alt exon, alt donor/acceptors, IR)
    td_data["flag_no_ovlp_nt_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["flag_no_ovlp_nt"] == 1),
        1,
        0,
    )
    td_data["flag_alt_exon_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["flag_alt_exon"] == 1), 1, 0
    )
    td_data["flag_alt_DA_RMP"] = np.where(
        (td_data["flag_RMP"] == 1)
        & (td_data["flag_alt_DA"] == 1),
        1,
        0,
    )
    td_data["flag_IR_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["flag_IR"] == 1), 1, 0
    )
    td_data["flag_5_var_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["flag_5_var"] == 1),
        1,
        0,
    )
    td_data["flag_3_var_RMP"] = np.where(
        (td_data["flag_RMP"] == 1) & (td_data["flag_3_var"] == 1),
        1,
        0,
    )
    # Get reciprocal minimum match gene category
    conditionsRecip = [
        (td_data["transcript_inGene"] == name1 + "_greater")
        & (td_data["num_RMP_inGene"] == 0),
        (td_data["transcript_inGene"] == name1 + "_greater")
        & (td_data["num_RMP_inGene"] > 0)
        & (
            td_data["num_RMP_inGene"]
            < td_data["num_transcript_inGene_" + name2]
        ),
        (td_data["transcript_inGene"] == name1 + "_greater")
        & (td_data["num_RMP_inGene"] > 0)
        & (
            td_data["num_RMP_inGene"]
            == td_data["num_transcript_inGene_" + name2]
        ),
        (td_data["transcript_inGene"] == name2 + "_greater")
        & (td_data["num_RMP_inGene"] == 0),
        (td_data["transcript_inGene"] == name2 + "_greater")
        & (td_data["num_RMP_inGene"] > 0)
        & (
            td_data["num_RMP_inGene"]
            < td_data["num_transcript_inGene_" + name1]
        ),
        (td_data["transcript_inGene"] == name2 + "_greater")
        & (td_data["num_RMP_inGene"] > 0)
        & (
            td_data["num_RMP_inGene"]
            == td_data["num_transcript_inGene_" + name1]
        ),
        (td_data["transcript_inGene"] == "match")
        & (td_data["num_RMP_inGene"] == 0),
        (td_data["transcript_inGene"] == "match")
        & (td_data["num_RMP_inGene"] > 0)
        & (
            td_data["num_RMP_inGene"]
            < td_data["num_transcript_inGene_" + name1]
        ),
        (td_data["transcript_inGene"] == "match")
        & (td_data["num_RMP_inGene"] > 0)
        & (
            td_data["num_RMP_inGene"]
            == td_data["num_transcript_inGene_" + name1]
        ),
    ]
    choicesRecrip = [
        "no_reciprocal_pairs",
        "partial_reciprocal_pairs",
        "reciprocal_pairs",
        "no_reciprocal_pairs",
        "partial_reciprocal_pairs",
        "reciprocal_pairs",
        "no_reciprocal_pairs",
        "partial_reciprocal_pairs",
        "reciprocal_pairs",
    ]
    td_data["recip_min_pair_inGene"] = np.select(
        conditionsRecip, choicesRecrip, "missing"
    )
    if len(td_data[td_data["recip_min_pair_inGene"] == "missing"]) > 0:
        logger.error(
            "Unexpected variable assignment for gene reciprocal minimum pair comparison"
        )

    # Flag minimum pair with considering ties (all ties are flagged as minimum)
    # First make variable numbering each group of unique distance values and T1 (or T2)
    #   (this means all T2 with the same distance to T1 will get the same group number)
    # Then groupby T1 (or T2) only and find the minimum group number
    #   (transcripts were ordered by min distance, so the min group # is the group with the min distance)
    # Then flag all rows where the group number is equal to the minimum group number
    #   (this will flag all T2 that share the same minimum distance to T1)
    # Remove intermediate columns and keep just the flags
    td_data["num_group_" + name1] = td_data.groupby(
        [
            "prop_ER_noOvlp",
            "prop_jxn_noOvlp",
            "prop_nt_noOvlp",
            "transcript_1"
        ])["transcript_2"].ngroup()
    td_data["min_group_num_" + name1] = td_data.groupby("transcript_1")[
            "num_group_" + name1].transform('min')
    td_data["flag_min_match_" + name1 + "_w_tie"] = np.where(
            td_data["min_group_num_" + name1] == td_data["num_group_" + name1],
            1,
            0
    )
    td_data["num_group_" + name2] = td_data.groupby(
        [
            "prop_ER_noOvlp",
            "prop_jxn_noOvlp",
            "prop_nt_noOvlp",
            "transcript_2"
        ])["transcript_2"].ngroup()
    td_data["min_group_num_" + name2] = td_data.groupby("transcript_2")[
            "num_group_" + name2].transform('min')
    td_data["flag_min_match_" + name2 + "_w_tie"] = np.where(
            td_data["min_group_num_" + name2] == td_data["num_group_" + name2],
            1,
            0
    )
    # Flag reciprocal minimum pairs with considering ties
    td_data["flag_RMP_w_tie"] = np.where(
        (td_data["flag_min_match_" + name1 + "_w_tie"] == 1)
        & (td_data["flag_min_match_" + name2 + "_w_tie"] == 1),
        1,
        0,
    )

    # !!! Check min distance with ties flags
    # check1/check2: Check that each T1 (or T2) has at least one minimum
    # check3/check4: Check that each T1 (or T2) that is flagged for a tie has more than one minimum
    # check5/check6: Check that for minimum (w/ ties) of each T1 (or T2) has the same number of
    #   flags as there are transcripts in the tie list
    #   (will be 1 if there was no tie)
    td_data["num_min_match_" + name1 + "_w_tie"] = td_data.groupby("transcript_1")[
            "flag_min_match_" + name1 + "_w_tie"].transform('sum')
    td_data["num_min_match_" + name2 + "_w_tie"] = td_data.groupby("transcript_2")[
            "flag_min_match_" + name2 + "_w_tie"].transform('sum')
    td_data["num_tie_" + name1] = td_data[
            name1 + "_distance_tie"].str.split("|").apply(lambda x: len(x))
    td_data["num_tie_" + name2] = td_data[
            name2 + "_distance_tie"].str.split("|").apply(lambda x: len(x))
    check1 = (td_data["num_min_match_" + name1 + "_w_tie"] > 0).all()
    check2 = (td_data["num_min_match_" + name2 + "_w_tie"] > 0).all()
    check3 = (td_data[
            (td_data["flag_" + name1 + "_tie"]==1)
            & (td_data["flag_min_match_" + name1 + "_w_tie"]==1)
        ]["num_min_match_" + name1 + "_w_tie"] > 1).all()
    check4 = (td_data[
            (td_data["flag_" + name2 + "_tie"]==1)
            & (td_data["flag_min_match_" + name2 + "_w_tie"]==1)
        ]["num_min_match_" + name2 + "_w_tie"] > 1).all()
    check5 = td_data[td_data["flag_min_match_" + name1 + "_w_tie"]==1
                     ]["num_min_match_" + name1 + "_w_tie"].equals(
                     td_data[td_data["flag_min_match_" + name1 + "_w_tie"]==1
                             ]["num_tie_" + name1])
    check6 = td_data[td_data["flag_min_match_" + name2 + "_w_tie"]==1
                     ]["num_min_match_" + name2 + "_w_tie"].equals(
                     td_data[td_data["flag_min_match_" + name2 + "_w_tie"]==1
                             ]["num_tie_" + name2])
    if not check1 or not check2 or not check3 or not check4 or not check5 or not check6:
        logger.error("Unexpected association of minimum distance within ties.")

    # Count number of reciprocal minimum pairs in each gene for T1 or T2
    # NOTE: These values could be different if one GTF has more transcripts than
    #   the other and the extra transcript matches another, leading to a tie.
    #   If both GTF have the same # of total transcripts and this equals the
    #   number of recip min for both T1 and T2 then this is a "match reciprocal minimum"
    #   gene. The maximal number of reciprocal minimum pairs is the minimum
    #   of the reciprocal minimums of T1 or T2
    # Flag T1 (or T2) if it has a reciprocal minimum with ties
    td_data["num_RMP_inGene_w_tie"] = get_num_inGene_w_tie(
            td_data,
            "flag_RMP_w_tie"
    )

    # For reciprocal minimums (with considering ties)
    # Flag FSM (all junctions matching) and ERS (all exonic regions shared) min matches
    td_data["flag_identical_RMP_w_tie"] = np.where(
        (td_data["flag_RMP_w_tie"] == 1) & (td_data["prop_nt_noOvlp"] == 0), 1, 0
    )
    td_data["flag_FSM_RMP_w_tie"] = np.where(
        (td_data["flag_RMP_w_tie"] == 1) & (td_data["prop_jxn_noOvlp"] == 0),
        1,
        0,
    )
    td_data["flag_ERS_RMP_w_tie"] = np.where(
        (td_data["flag_RMP_w_tie"] == 1) & (td_data["prop_ER_noOvlp"] == 0), 1, 0
    )
    td_data["flag_ERS_noIR_RMP_w_tie"] = np.where(
        (td_data["flag_RMP_w_tie"] == 1)
        & (td_data["prop_ER_noOvlp"] == 0)
        & (td_data["num_IR_EF_T1"] + td_data["num_IR_EF_T2"] == 0),
        1,
        0,
    )
    td_data["flag_ERS_wIR_RMP_w_tie"] = np.where(
        (td_data["flag_RMP_w_tie"] == 1)
        & (td_data["prop_ER_noOvlp"] == 0)
        & (td_data["num_IR_EF_T1"] + td_data["num_IR_EF_T2"] > 0),
        1,
        0,
    )
    td_data["num_identical_RMP_inGene_w_tie"] = get_num_inGene_w_tie(
            td_data,
            "flag_identical_RMP_w_tie"
    )
    td_data["num_FSM_RMP_inGene_w_tie"] = get_num_inGene_w_tie(
            td_data,
            "flag_FSM_RMP_w_tie"
    )
    td_data["num_ERS_RMP_inGene_w_tie"] = get_num_inGene_w_tie(
            td_data,
            "flag_ERS_RMP_w_tie"
    )
    td_data["num_ERS_noIR_RMP_inGene_w_tie"] = get_num_inGene_w_tie(
            td_data,
            "flag_ERS_noIR_RMP_w_tie"
    )
    td_data["num_ERS_wIR_RMP_inGene_w_tie"] = get_num_inGene_w_tie(
            td_data,
            "flag_ERS_wIR_RMP_w_tie"
    )
    td_data["prop_identical_RMP_inGene_w_tie"] = (
        td_data["num_identical_RMP_inGene_w_tie"]
        / td_data["num_RMP_inGene_w_tie"]
    )
    td_data["prop_FSM_RMP_inGene_w_tie"] = (
        td_data["num_FSM_RMP_inGene_w_tie"]
        / td_data["num_RMP_inGene_w_tie"]
    )
    td_data["prop_ERS_RMP_inGene_w_tie"] = (
        td_data["num_ERS_RMP_inGene_w_tie"]
        / td_data["num_RMP_inGene_w_tie"]
    )
    td_data["prop_ERS_noIR_RMP_inGene_w_tie"] = (
        td_data["num_ERS_noIR_RMP_inGene_w_tie"]
        / td_data["num_RMP_inGene_w_tie"]
    )
    td_data["prop_ERS_wIR_RMP_inGene_w_tie"] = (
        td_data["num_ERS_wIR_RMP_inGene_w_tie"]
        / td_data["num_RMP_inGene_w_tie"]
    )

    # Get reciprocal minimum match gene category
    conditionsRecipTies = [
        (td_data["transcript_inGene"] == name1 + "_greater")
        & (td_data["num_RMP_inGene_w_tie"] == 0),
        (td_data["transcript_inGene"] == name1 + "_greater")
        & (td_data["num_RMP_inGene_w_tie"] > 0)
        & (
            td_data["num_RMP_inGene_w_tie"]
            < td_data["num_transcript_inGene_" + name2]
        ),
        (td_data["transcript_inGene"] == name1 + "_greater")
        & (td_data["num_RMP_inGene_w_tie"] > 0)
        & (
            td_data["num_RMP_inGene_w_tie"]
            == td_data["num_transcript_inGene_" + name2]
        ),
        (td_data["transcript_inGene"] == name2 + "_greater")
        & (td_data["num_RMP_inGene_w_tie"] == 0),
        (td_data["transcript_inGene"] == name2 + "_greater")
        & (td_data["num_RMP_inGene_w_tie"] > 0)
        & (
            td_data["num_RMP_inGene_w_tie"]
            < td_data["num_transcript_inGene_" + name1]
        ),
        (td_data["transcript_inGene"] == name2 + "_greater")
        & (td_data["num_RMP_inGene_w_tie"] > 0)
        & (
            td_data["num_RMP_inGene_w_tie"]
            == td_data["num_transcript_inGene_" + name1]
        ),
        (td_data["transcript_inGene"] == "match")
        & (td_data["num_RMP_inGene_w_tie"] == 0),
        (td_data["transcript_inGene"] == "match")
        & (td_data["num_RMP_inGene_w_tie"] > 0)
        & (
            td_data["num_RMP_inGene_w_tie"]
            < td_data["num_transcript_inGene_" + name1]
        ),
        (td_data["transcript_inGene"] == "match")
        & (td_data["num_RMP_inGene_w_tie"] > 0)
        & (
            td_data["num_RMP_inGene_w_tie"]
            == td_data["num_transcript_inGene_" + name1]
        ),
    ]
    choicesRecripTies = [
        "no_reciprocal_pairs",
        "partial_reciprocal_pairs",
        "reciprocal_pairs",
        "no_reciprocal_pairs",
        "partial_reciprocal_pairs",
        "reciprocal_pairs",
        "no_reciprocal_pairs",
        "partial_reciprocal_pairs",
        "reciprocal_pairs",
    ]
    td_data["recip_min_pair_inGene_w_tie"] = np.select(
        conditionsRecipTies, choicesRecripTies, "missing"
    )
    if len(td_data[td_data["recip_min_pair_inGene_w_tie"] == "missing"]) > 0:
        logger.error(
            "Unexpected variable assignment for gene reciprocal minimum pair comparison"
        )
        
    # Remove intermediate columns
    td_data = td_data.drop(columns=[
            "num_group_" + name1,
            "min_group_num_" + name1,
            "num_group_" + name2,
            "min_group_num_" + name2,
            "num_min_match_" + name1 + "_w_tie",
            "num_min_match_" + name2 + "_w_tie",
            "num_tie_" + name1,
            "num_tie_" + name2,
    ])

    # Return minimum distance of transcript all pairs
    if out_pairs == "all":
        return td_data
    # Return minimum distance of only minimum pairs of each transcript
    elif out_pairs == "both":
        return td_data[
            (td_data["flag_min_match_" + name1 + "_w_tie"] == 1)
            | (td_data["flag_min_match_" + name2 + "_w_tie"] == 1)
        ]
    elif out_pairs == "first":
        return td_data[td_data["flag_min_match_" + name1 + "_w_tie"] == 1]
    elif out_pairs == "second":
        return td_data[td_data["flag_min_match_" + name2 + "_w_tie"] == 1]
    else:
        logger.error("Unexpected value for output pairs argument.")

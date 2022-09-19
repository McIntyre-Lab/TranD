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
            'num_junction_T1_only',
            'num_junction_T2_only',
            'num_junction_shared',
            'prop_junction_diff',
            'prop_junction_similar',
            'junction_T1_all',
            'junction_T2_all',
            'junction_T1_only',
            'junction_T2_only',
            'junction_shared',
            'num_ER_T1_only',
            'num_ER_T2_only',
            'num_ER_shared',
            'prop_ER_diff',
            'prop_ER_similar',
            'ER_T1_only',
            'ER_T2_only',
            'ER_shared',
            'num_fragment_T1_only',
            'num_fragment_T2_only',
            'num_fragment_shared',
            'prop_fragment_diff',
            'prop_fragment_similar',
            'fragment_T1_only',
            'fragment_T2_only',
            'fragment_shared',
            'num_fragment_singletons_T1_only',
            'num_fragment_singletons_T2_only',
            'num_fragment_singletons_shared',
            'num_IR_fragment_T1',
            'num_IR_fragment_T2',
            'IR_fragment_T1',
            'IR_fragment_T2',
            'num_nt_shared',
            'num_nt_T1_only',
            'num_nt_T2_only',
            'num_nt_diff',
            'total_nt',
            'prop_nt_diff',
            'prop_nt_similar',
            'num_nt_T1_only_in_shared_ER',
            'num_nt_T2_only_in_shared_ER',
            'num_nt_shared_in_shared_ER',
            'total_nt_in_shared_ER',
            'prop_nt_diff_in_shared_ER',
            'prop_nt_similar_in_shared_ER',
            'num_nt_T1_only_in_unique_ER',
            'num_nt_T2_only_in_unique_ER',
            'flag_FSM',
            'flag_T1_ISM_of_T2',
            'flag_T2_ISM_of_T1',
            'flag_IR',
            'flag_5_variation',
            'flag_3_variation',
            'flag_alt_donor_acceptor',
            'flag_alt_exon',
            'flag_no_shared_nt',
            'num_transcript_in_gene_'+name1,
            'num_transcript_in_gene_'+name2,
            'flag_' + name1 + '_greater',
            'flag_' + name2 + '_greater',
            'flag_match',
            'transcript_in_gene',
            'min_match_' + name1,
            'flag_min_match_' + name1,
            'min_match_' + name2,
            'flag_min_match_' + name2,
            'flag_recip_min_match',
            'num_recip_min_match_in_gene',
            'flag_' + name1 + '_tie',
            'flag_' + name2 + '_tie',
            name1 + '_distance_ties',
            name2 + '_distance_ties',
            'flag_identical_recip_min_match',
            'flag_FSM_recip_min_match',
            'flag_ERM_recip_min_match',
            'flag_ERM_noIR_recip_min_match',
            'flag_ERM_withIR_recip_min_match',
            'num_identical_recip_min_match_in_gene',
            'num_FSM_recip_min_match_in_gene',
            'num_ERM_recip_min_match_in_gene',
            'num_ERM_noIR_recip_min_match_in_gene',
            'num_ERM_withIR_recip_min_match_in_gene',
            'prop_identical_recip_min_match_in_gene',
            'prop_FSM_recip_min_match_in_gene',
            'prop_ERM_recip_min_match_in_gene',
            'prop_ERM_noIR_recip_min_match_in_gene',
            'prop_ERM_withIR_recip_min_match_in_gene',
            'flag_no_shared_nt_recip_min_match',
            'flag_alt_exon_recip_min_match',
            'flag_alt_donor_acceptor_recip_min_match',
            'flag_IR_recip_min_match',
            'flag_5_variation_recip_min_match',
            'flag_3_variation_recip_min_match',
            'recip_min_pair_in_gene',
            'flag_min_match_' + name1 + '_w_ties',
            'flag_min_match_' + name2 + '_w_ties',
            'flag_recip_min_match_w_ties',
            'num_recip_min_match_in_gene_w_ties',
            'flag_identical_recip_min_match_w_ties',
            'flag_FSM_recip_min_match_w_ties',
            'flag_ERM_recip_min_match_w_ties',
            'flag_ERM_noIR_recip_min_match_w_ties',
            'flag_ERM_withIR_recip_min_match_w_ties',
            'num_identical_recip_min_match_in_gene_w_ties',
            'num_FSM_recip_min_match_in_gene_w_ties',
            'num_ERM_recip_min_match_in_gene_w_ties',
            'num_ERM_noIR_recip_min_match_in_gene_w_ties',
            'prop_identical_recip_min_match_in_gene_w_ties',
            'prop_FSM_recip_min_match_in_gene_w_ties',
            'prop_ERM_recip_min_match_in_gene_w_ties',
            'prop_ERM_noIR_recip_min_match_in_gene_w_ties',
            'num_ERM_withIR_recip_min_match_in_gene_w_ties',
            'prop_ERM_withIR_recip_min_match_in_gene_w_ties',
            'recip_min_pair_in_gene_w_ties',
        ]
    return md_df_cols

def get_num_in_gene_w_ties(df, flag_name):
    df2 = df.copy()
    df2["name_1_flag_w_ties"] = np.where(
           df2[flag_name] == 1,
           df2["transcript_1"],
           np.nan
    )
    df2["name_2_flag_w_ties"] = np.where(
           df2[flag_name] == 1,
           df2["transcript_2"],
           np.nan
    )
    df2["num_name_1_flag_w_ties"] = (
        df2.groupby("gene_id")[
            "name_1_flag_w_ties"
        ].transform("nunique")
    )
    df2["num_name_2_flag_w_ties"] = (
        df2.groupby("gene_id")[
            "name_2_flag_w_ties"
        ].transform("nunique")
    )
    df2["num_identical_recip_min_match_in_gene_w_ties"] = df2[[
            "num_name_1_flag_w_ties",
            "num_name_2_flag_w_ties"
        ]].min(axis=1)
    return df2["num_identical_recip_min_match_in_gene_w_ties"]

def identify_min_pair(td_data, out_pairs, name1, name2):

    # Get number of transcripts per gene
    td_data["num_transcript_in_gene_" + name1] = td_data.groupby("gene_id")[
        "transcript_1"
    ].transform("nunique")
    td_data["num_transcript_in_gene_" + name2] = td_data.groupby("gene_id")[
        "transcript_2"
    ].transform("nunique")

    # Flag groups based on transcripts per gene
    td_data["flag_" + name1 + "_greater"] = np.where(
        (
            td_data["num_transcript_in_gene_" + name1]
            > td_data["num_transcript_in_gene_" + name2]
        )
        & (td_data["num_transcript_in_gene_" + name1] > 0)
        & (td_data["num_transcript_in_gene_" + name2] > 0),
        1,
        0,
    )
    td_data["flag_" + name2 + "_greater"] = np.where(
        (
            td_data["num_transcript_in_gene_" + name1]
            < td_data["num_transcript_in_gene_" + name2]
        )
        & (td_data["num_transcript_in_gene_" + name1] > 0)
        & (td_data["num_transcript_in_gene_" + name2] > 0),
        1,
        0,
    )
    td_data["flag_match"] = np.where(
        (
            td_data["num_transcript_in_gene_" + name1]
            == td_data["num_transcript_in_gene_" + name2]
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
    td_data["transcript_in_gene"] = np.select(conditionsGene, choicesGene, "missing")
    if len(td_data[td_data["transcript_in_gene"] == "missing"]) > 0:
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
            "prop_ER_diff",
            "prop_junction_diff",
            "prop_nt_diff",
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
    td_data["flag_recip_min_match"] = np.where(
        (td_data["flag_min_match_" + name1] == 1)
        & (td_data["flag_min_match_" + name2] == 1),
        1,
        0,
    )
    td_data["num_recip_min_match_in_gene"] = td_data.groupby("gene_id")[
        "flag_recip_min_match"
    ].transform("sum")

    # Flag pairs that are ties in between at least 2 transcirpts in the other dataset
    td_data["flag_" + name1 + "_tie"] = (
        td_data[
            [
                "transcript_1",
                "prop_junction_diff",
                "prop_ER_diff",
                "prop_nt_diff",
                "num_nt_diff",
            ]
        ]
        .duplicated(keep=False)
        .astype(int)
    )
    td_data["flag_" + name2 + "_tie"] = (
        td_data[
            [
                "transcript_2",
                "prop_junction_diff",
                "prop_ER_diff",
                "prop_nt_diff",
                "num_nt_diff",
            ]
        ]
        .duplicated(keep=False)
        .astype(int)
    )

    # List the other transcripts that the distances are tied with
    td_data[name1 + "_distance_ties"] = td_data.groupby(
            [
                "transcript_1",
                "prop_junction_diff",
                "prop_ER_diff",
                "prop_nt_diff",
                "num_nt_diff",
            ])["transcript_2"].transform(lambda x: "|".join(x))
    td_data[name2 + "_distance_ties"] = td_data.groupby(
            [
                "transcript_2",
                "prop_junction_diff",
                "prop_ER_diff",
                "prop_nt_diff",
                "num_nt_diff",
            ])["transcript_1"].transform(lambda x: "|".join(x))
        
    # Flag FSM (all junctions matching) and ERM (all exonic regions shared) min matches
    td_data["flag_identical_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["prop_nt_diff"] == 0),
        1,
        0
    )
    td_data["flag_FSM_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["prop_junction_diff"] == 0),
        1,
        0,
    )
    td_data["flag_ERM_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["prop_ER_diff"] == 0),
        1,
        0
    )
    td_data["flag_ERM_noIR_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1)
        & (td_data["prop_ER_diff"] == 0)
        & (td_data["num_IR_fragment_T1"] + td_data["num_IR_fragment_T2"] == 0),
        1,
        0,
    )
    td_data["flag_ERM_withIR_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1)
        & (td_data["prop_ER_diff"] == 0)
        & (td_data["num_IR_fragment_T1"] + td_data["num_IR_fragment_T2"] > 0),
        1,
        0,
    )
    td_data["num_identical_recip_min_match_in_gene"] = td_data.groupby("gene_id")[
        "flag_identical_recip_min_match"
    ].transform("sum")
    td_data["num_FSM_recip_min_match_in_gene"] = td_data.groupby("gene_id")[
        "flag_FSM_recip_min_match"
    ].transform("sum")
    td_data["num_ERM_recip_min_match_in_gene"] = td_data.groupby("gene_id")[
        "flag_ERM_recip_min_match"
    ].transform("sum")
    td_data["num_ERM_noIR_recip_min_match_in_gene"] = td_data.groupby("gene_id")[
        "flag_ERM_noIR_recip_min_match"
    ].transform("sum")
    td_data["num_ERM_withIR_recip_min_match_in_gene"] = td_data.groupby("gene_id")[
        "flag_ERM_withIR_recip_min_match"
    ].transform("sum")
    td_data["prop_identical_recip_min_match_in_gene"] = (
        td_data["num_identical_recip_min_match_in_gene"]
        / td_data["num_recip_min_match_in_gene"]
    )
    td_data["prop_FSM_recip_min_match_in_gene"] = (
        td_data["num_FSM_recip_min_match_in_gene"]
        / td_data["num_recip_min_match_in_gene"]
    )
    td_data["prop_ERM_recip_min_match_in_gene"] = (
        td_data["num_ERM_recip_min_match_in_gene"]
        / td_data["num_recip_min_match_in_gene"]
    )
    td_data["prop_ERM_noIR_recip_min_match_in_gene"] = (
        td_data["num_ERM_noIR_recip_min_match_in_gene"]
        / td_data["num_recip_min_match_in_gene"]
    )
    td_data["prop_ERM_withIR_recip_min_match_in_gene"] = (
        td_data["num_ERM_withIR_recip_min_match_in_gene"]
        / td_data["num_recip_min_match_in_gene"]
    )

    # Flag different types of splicing differences in min matches
    #   (alt exon, alt donor/acceptors, IR)
    td_data["flag_no_shared_nt_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["flag_no_shared_nt"] == 1),
        1,
        0,
    )
    td_data["flag_alt_exon_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["flag_alt_exon"] == 1), 1, 0
    )
    td_data["flag_alt_donor_acceptor_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1)
        & (td_data["flag_alt_donor_acceptor"] == 1),
        1,
        0,
    )
    td_data["flag_IR_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["flag_IR"] == 1), 1, 0
    )
    td_data["flag_5_variation_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["flag_5_variation"] == 1),
        1,
        0,
    )
    td_data["flag_3_variation_recip_min_match"] = np.where(
        (td_data["flag_recip_min_match"] == 1) & (td_data["flag_3_variation"] == 1),
        1,
        0,
    )
    # Get reciprocal minimum match gene category
    conditionsRecip = [
        (td_data["transcript_in_gene"] == name1 + "_greater")
        & (td_data["num_recip_min_match_in_gene"] == 0),
        (td_data["transcript_in_gene"] == name1 + "_greater")
        & (td_data["num_recip_min_match_in_gene"] > 0)
        & (
            td_data["num_recip_min_match_in_gene"]
            < td_data["num_transcript_in_gene_" + name2]
        ),
        (td_data["transcript_in_gene"] == name1 + "_greater")
        & (td_data["num_recip_min_match_in_gene"] > 0)
        & (
            td_data["num_recip_min_match_in_gene"]
            == td_data["num_transcript_in_gene_" + name2]
        ),
        (td_data["transcript_in_gene"] == name2 + "_greater")
        & (td_data["num_recip_min_match_in_gene"] == 0),
        (td_data["transcript_in_gene"] == name2 + "_greater")
        & (td_data["num_recip_min_match_in_gene"] > 0)
        & (
            td_data["num_recip_min_match_in_gene"]
            < td_data["num_transcript_in_gene_" + name1]
        ),
        (td_data["transcript_in_gene"] == name2 + "_greater")
        & (td_data["num_recip_min_match_in_gene"] > 0)
        & (
            td_data["num_recip_min_match_in_gene"]
            == td_data["num_transcript_in_gene_" + name1]
        ),
        (td_data["transcript_in_gene"] == "match")
        & (td_data["num_recip_min_match_in_gene"] == 0),
        (td_data["transcript_in_gene"] == "match")
        & (td_data["num_recip_min_match_in_gene"] > 0)
        & (
            td_data["num_recip_min_match_in_gene"]
            < td_data["num_transcript_in_gene_" + name1]
        ),
        (td_data["transcript_in_gene"] == "match")
        & (td_data["num_recip_min_match_in_gene"] > 0)
        & (
            td_data["num_recip_min_match_in_gene"]
            == td_data["num_transcript_in_gene_" + name1]
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
    td_data["recip_min_pair_in_gene"] = np.select(
        conditionsRecip, choicesRecrip, "missing"
    )
    if len(td_data[td_data["recip_min_pair_in_gene"] == "missing"]) > 0:
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
            "prop_ER_diff",
            "prop_junction_diff",
            "prop_nt_diff",
            "transcript_1"
        ])["transcript_2"].ngroup()
    td_data["min_group_num_" + name1] = td_data.groupby("transcript_1")[
            "num_group_" + name1].transform('min')
    td_data["flag_min_match_" + name1 + "_w_ties"] = np.where(
            td_data["min_group_num_" + name1] == td_data["num_group_" + name1],
            1,
            0
    )
    td_data["num_group_" + name2] = td_data.groupby(
        [
            "prop_ER_diff",
            "prop_junction_diff",
            "prop_nt_diff",
            "transcript_2"
        ])["transcript_2"].ngroup()
    td_data["min_group_num_" + name2] = td_data.groupby("transcript_2")[
            "num_group_" + name2].transform('min')
    td_data["flag_min_match_" + name2 + "_w_ties"] = np.where(
            td_data["min_group_num_" + name2] == td_data["num_group_" + name2],
            1,
            0
    )
    # Flag reciprocal minimum pairs with considering ties
    td_data["flag_recip_min_match_w_ties"] = np.where(
        (td_data["flag_min_match_" + name1 + "_w_ties"] == 1)
        & (td_data["flag_min_match_" + name2 + "_w_ties"] == 1),
        1,
        0,
    )

    # !!! Check min distance with ties flags
    # check1/check2: Check that each T1 (or T2) has at least one minimum
    # check3/check4: Check that each T1 (or T2) that is flagged for a tie has more than one minimum
    # check5/check6: Check that for minimum (w/ ties) of each T1 (or T2) has the same number of
    #   flags as there are transcripts in the tie list
    #   (will be 1 if there was no tie)
    td_data["num_min_match_" + name1 + "_w_ties"] = td_data.groupby("transcript_1")[
            "flag_min_match_" + name1 + "_w_ties"].transform('sum')
    td_data["num_min_match_" + name2 + "_w_ties"] = td_data.groupby("transcript_2")[
            "flag_min_match_" + name2 + "_w_ties"].transform('sum')
    td_data["num_ties_" + name1] = td_data[
            name1 + "_distance_ties"].str.split("|").apply(lambda x: len(x))
    td_data["num_ties_" + name2] = td_data[
            name2 + "_distance_ties"].str.split("|").apply(lambda x: len(x))
    check1 = (td_data["num_min_match_" + name1 + "_w_ties"] > 0).all()
    check2 = (td_data["num_min_match_" + name2 + "_w_ties"] > 0).all()
    check3 = (td_data[
            (td_data["flag_" + name1 + "_tie"]==1)
            & (td_data["flag_min_match_" + name1 + "_w_ties"]==1)
        ]["num_min_match_" + name1 + "_w_ties"] > 1).all()
    check4 = (td_data[
            (td_data["flag_" + name2 + "_tie"]==1)
            & (td_data["flag_min_match_" + name2 + "_w_ties"]==1)
        ]["num_min_match_" + name2 + "_w_ties"] > 1).all()
    check5 = td_data[td_data["flag_min_match_" + name1 + "_w_ties"]==1
                     ]["num_min_match_" + name1 + "_w_ties"].equals(
                     td_data[td_data["flag_min_match_" + name1 + "_w_ties"]==1
                             ]["num_ties_" + name1])
    check6 = td_data[td_data["flag_min_match_" + name2 + "_w_ties"]==1
                     ]["num_min_match_" + name2 + "_w_ties"].equals(
                     td_data[td_data["flag_min_match_" + name2 + "_w_ties"]==1
                             ]["num_ties_" + name2])
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
    td_data["num_recip_min_match_in_gene_w_ties"] = get_num_in_gene_w_ties(
            td_data,
            "flag_recip_min_match_w_ties"
    )

    # For reciprocal minimums (with considering ties)
    # Flag FSM (all junctions matching) and ERM (all exonic regions shared) min matches
    td_data["flag_identical_recip_min_match_w_ties"] = np.where(
        (td_data["flag_recip_min_match_w_ties"] == 1) & (td_data["prop_nt_diff"] == 0), 1, 0
    )
    td_data["flag_FSM_recip_min_match_w_ties"] = np.where(
        (td_data["flag_recip_min_match_w_ties"] == 1) & (td_data["prop_junction_diff"] == 0),
        1,
        0,
    )
    td_data["flag_ERM_recip_min_match_w_ties"] = np.where(
        (td_data["flag_recip_min_match_w_ties"] == 1) & (td_data["prop_ER_diff"] == 0), 1, 0
    )
    td_data["flag_ERM_noIR_recip_min_match_w_ties"] = np.where(
        (td_data["flag_recip_min_match_w_ties"] == 1)
        & (td_data["prop_ER_diff"] == 0)
        & (td_data["num_IR_fragment_T1"] + td_data["num_IR_fragment_T2"] == 0),
        1,
        0,
    )
    td_data["flag_ERM_withIR_recip_min_match_w_ties"] = np.where(
        (td_data["flag_recip_min_match_w_ties"] == 1)
        & (td_data["prop_ER_diff"] == 0)
        & (td_data["num_IR_fragment_T1"] + td_data["num_IR_fragment_T2"] > 0),
        1,
        0,
    )
    td_data["num_identical_recip_min_match_in_gene_w_ties"] = get_num_in_gene_w_ties(
            td_data,
            "flag_identical_recip_min_match_w_ties"
    )
    td_data["num_FSM_recip_min_match_in_gene_w_ties"] = get_num_in_gene_w_ties(
            td_data,
            "flag_FSM_recip_min_match_w_ties"
    )
    td_data["num_ERM_recip_min_match_in_gene_w_ties"] = get_num_in_gene_w_ties(
            td_data,
            "flag_ERM_recip_min_match_w_ties"
    )
    td_data["num_ERM_noIR_recip_min_match_in_gene_w_ties"] = get_num_in_gene_w_ties(
            td_data,
            "flag_ERM_noIR_recip_min_match_w_ties"
    )
    td_data["num_ERM_withIR_recip_min_match_in_gene_w_ties"] = get_num_in_gene_w_ties(
            td_data,
            "flag_ERM_withIR_recip_min_match_w_ties"
    )
    td_data["prop_identical_recip_min_match_in_gene_w_ties"] = (
        td_data["num_identical_recip_min_match_in_gene_w_ties"]
        / td_data["num_recip_min_match_in_gene_w_ties"]
    )
    td_data["prop_FSM_recip_min_match_in_gene_w_ties"] = (
        td_data["num_FSM_recip_min_match_in_gene_w_ties"]
        / td_data["num_recip_min_match_in_gene_w_ties"]
    )
    td_data["prop_ERM_recip_min_match_in_gene_w_ties"] = (
        td_data["num_ERM_recip_min_match_in_gene_w_ties"]
        / td_data["num_recip_min_match_in_gene_w_ties"]
    )
    td_data["prop_ERM_noIR_recip_min_match_in_gene_w_ties"] = (
        td_data["num_ERM_noIR_recip_min_match_in_gene_w_ties"]
        / td_data["num_recip_min_match_in_gene_w_ties"]
    )
    td_data["prop_ERM_withIR_recip_min_match_in_gene_w_ties"] = (
        td_data["num_ERM_withIR_recip_min_match_in_gene_w_ties"]
        / td_data["num_recip_min_match_in_gene_w_ties"]
    )

    # Get reciprocal minimum match gene category
    conditionsRecipTies = [
        (td_data["transcript_in_gene"] == name1 + "_greater")
        & (td_data["num_recip_min_match_in_gene_w_ties"] == 0),
        (td_data["transcript_in_gene"] == name1 + "_greater")
        & (td_data["num_recip_min_match_in_gene_w_ties"] > 0)
        & (
            td_data["num_recip_min_match_in_gene_w_ties"]
            < td_data["num_transcript_in_gene_" + name2]
        ),
        (td_data["transcript_in_gene"] == name1 + "_greater")
        & (td_data["num_recip_min_match_in_gene_w_ties"] > 0)
        & (
            td_data["num_recip_min_match_in_gene_w_ties"]
            == td_data["num_transcript_in_gene_" + name2]
        ),
        (td_data["transcript_in_gene"] == name2 + "_greater")
        & (td_data["num_recip_min_match_in_gene_w_ties"] == 0),
        (td_data["transcript_in_gene"] == name2 + "_greater")
        & (td_data["num_recip_min_match_in_gene_w_ties"] > 0)
        & (
            td_data["num_recip_min_match_in_gene_w_ties"]
            < td_data["num_transcript_in_gene_" + name1]
        ),
        (td_data["transcript_in_gene"] == name2 + "_greater")
        & (td_data["num_recip_min_match_in_gene_w_ties"] > 0)
        & (
            td_data["num_recip_min_match_in_gene_w_ties"]
            == td_data["num_transcript_in_gene_" + name1]
        ),
        (td_data["transcript_in_gene"] == "match")
        & (td_data["num_recip_min_match_in_gene_w_ties"] == 0),
        (td_data["transcript_in_gene"] == "match")
        & (td_data["num_recip_min_match_in_gene_w_ties"] > 0)
        & (
            td_data["num_recip_min_match_in_gene_w_ties"]
            < td_data["num_transcript_in_gene_" + name1]
        ),
        (td_data["transcript_in_gene"] == "match")
        & (td_data["num_recip_min_match_in_gene_w_ties"] > 0)
        & (
            td_data["num_recip_min_match_in_gene_w_ties"]
            == td_data["num_transcript_in_gene_" + name1]
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
    td_data["recip_min_pair_in_gene_w_ties"] = np.select(
        conditionsRecipTies, choicesRecripTies, "missing"
    )
    if len(td_data[td_data["recip_min_pair_in_gene_w_ties"] == "missing"]) > 0:
        logger.error(
            "Unexpected variable assignment for gene reciprocal minimum pair comparison"
        )
        
    # Remove intermediate columns
    td_data = td_data.drop(columns=[
            "num_group_" + name1,
            "min_group_num_" + name1,
            "num_group_" + name2,
            "min_group_num_" + name2,
            "num_min_match_" + name1 + "_w_ties",
            "num_min_match_" + name2 + "_w_ties",
            "num_ties_" + name1,
            "num_ties_" + name2,
    ])

    # Return minimum distance of transcript all pairs
    if out_pairs == "all":
        return td_data
    # Return minimum distance of only minimum pairs of each transcript
    elif out_pairs == "both":
        return td_data[
            (td_data["flag_min_match_" + name1 + "_w_ties"] == 1)
            | (td_data["flag_min_match_" + name2 + "_w_ties"] == 1)
        ]
    elif out_pairs == "first":
        return td_data[td_data["flag_min_match_" + name1 + "_w_ties"] == 1]
    elif out_pairs == "second":
        return td_data[td_data["flag_min_match_" + name2 + "_w_ties"] == 1]
    else:
        logger.error("Unexpected value for output pairs argument.")

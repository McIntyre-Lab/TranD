#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:03:09 2021

@author: adalena.nanni
"""

import pandas as pd
import numpy as np


td_df_cols = [
    "gene_id",
    "transcript_1",
    "transcript_2",
    "num_junction_T1_only",
    "num_junction_T2_only",
    "num_junction_shared",
    "prop_junction_diff",
    "prop_junction_similar",
    "junction_T1_all",
    "junction_T2_all",
    "junction_T1_only",
    "junction_T2_only",
    "junction_shared",
    "num_ER_T1_only",
    "num_ER_T2_only",
    "num_ER_shared",
    "prop_ER_diff",
    "prop_ER_similar",
    "ER_T1_only",
    "ER_T2_only",
    "ER_shared",
    "num_fragment_T1_only",
    "num_fragment_T2_only",
    "num_fragment_shared",
    "prop_fragment_diff",
    "prop_fragment_similar",
    "fragment_T1_only",
    "fragment_T2_only",
    "fragment_shared",
    "num_fragment_singletons_T1_only",
    "num_fragment_singletons_T2_only",
    "num_fragment_singletons_shared",
    "num_IR_fragment_T1",
    "num_IR_fragment_T2",
    "IR_fragment_T1",
    "IR_fragment_T2",
    "num_nt_shared",
    "num_nt_T1_only",
    "num_nt_T2_only",
    "num_nt_diff",
    "total_nt",
    "prop_nt_diff",
    "prop_nt_similar",
    "num_nt_T1_only_in_shared_ER",
    "num_nt_T2_only_in_shared_ER",
    "num_nt_shared_in_shared_ER",
    "total_nt_in_shared_ER",
    "prop_nt_diff_in_shared_ER",
    "prop_nt_similar_in_shared_ER",
    "num_nt_T1_only_in_unique_ER",
    "num_nt_T2_only_in_unique_ER",
    "flag_FSM",
    "flag_T1_ISM_of_T2",
    "flag_T2_ISM_of_T1",
    "flag_IR",
    "flag_5_variation",
    "flag_3_variation",
    "flag_alt_donor_acceptor",
    "flag_alt_exon",
    "flag_no_shared_nt",
]


def calculate_distance(out_df, junction_df, gene_id, tx1_name, tx2_name, fsm=False):
    # Check for monoexons
    has_junctions = True
    if junction_df.empty:
        has_junctions = False
    # Set EF variable
    ef_df = out_df.copy()

    # Make empty Series to place distance values
    singlePair = pd.Series(index=td_df_cols, dtype=object)
    singlePair[["gene_id", "transcript_1", "transcript_2"]] = (
        gene_id,
        tx1_name,
        tx2_name,
    )

    # Set exon fragment and exon region ID values based on genomic coordinates
    ef_df["fragment_id"] = (
        ef_df["ef_chr"].map(str)
        + ":"
        + ef_df["ef_start"].map(str)
        + ":"
        + ef_df["ef_end"].map(str)
        + ":"
        + ef_df["ef_strand"].map(str)
    )
    ef_df["region_id"] = (
        ef_df["er_chr"].map(str)
        + ":"
        + ef_df["er_start"].map(str)
        + ":"
        + ef_df["er_end"].map(str)
        + ":"
        + ef_df["er_strand"].map(str)
    )

    # Flag shared fragments
    ef_df["flag_shared"] = np.where(
        (ef_df["transcript_id"] == tx1_name + "|" + tx2_name)
        | (ef_df["transcript_id"] == tx2_name + "|" + tx1_name),
        1,
        0,
    )

    # Flag singleton fragments (where the fragment is a full exon region)
    ef_df["num_fragment_in_ER"] = ef_df.groupby("er_id")["ef_id"].transform("count")
    ef_df["flag_singleton"] = np.where(ef_df["num_fragment_in_ER"] == 1, 1, 0)

    # Add fragment lengths (end - start)
    ef_df["fragment_length"] = ef_df["ef_end"].map(int) - ef_df["ef_start"].map(int)

    # Sort junctions by chromosome and coordinates
    if has_junctions:
        sorted_junction_df = junction_df.copy()
        sorted_junction_df[["chr", "start", "end", "strand"]] = (
                sorted_junction_df["coords"].str.split(":", expand=True)
        )
        sorted_junction_df[["start", "end"]] = (
                sorted_junction_df[["start", "end"]].astype(int)
        )
        sorted_junction_df = sorted_junction_df.sort_values(
                by=["chr", "start", "end"]
            ).reset_index(drop=True)
    else:
        sorted_junction_df = pd.DataFrame()

    # Get distance measures for junctions
    singlePair = get_junction_distance(
        singlePair, sorted_junction_df, tx1_name, tx2_name, fsm=fsm
    )

    # Gt distance measures for exon regions (ER) and exon fragments (EF)
    singlePair, ERSharedSet = get_ER_distance(
        singlePair, ef_df, tx1_name, tx2_name, fsm=fsm
    )
    singlePair = get_EF_distance(singlePair, ef_df, tx1_name, tx2_name, ERSharedSet)

    # Set flags for different alternative splicing (AS) events
    singlePair = set_AS_flags(singlePair, ef_df, tx1_name, tx2_name, has_junctions)

    # Return distance of transcript pair
    return singlePair


def get_junction_distance(singlePair, sorted_junction_df, tx1_name, tx2_name, fsm=False):
    # Check if transcript pair shares all junctions (FSM, full-splice match)
    if fsm:
        # Check if both transcripts are monoexon (no junctions)
        if len(sorted_junction_df) == 0:
            singlePair[
                ["num_junction_T1_only", "num_junction_T2_only", "num_junction_shared"]
            ] = 0
            singlePair["prop_junction_diff"] = 0
            singlePair["prop_junction_similar"] = 1
            singlePair[
                ["junction_T1_only",
                 "junction_T2_only",
                 "junction_shared",
                 "junction_T1_all",
                 "junction_T2_all"]
            ] = ""
        # FSM transcripts are multiexon
        else:
            singlePair[["num_junction_T1_only", "num_junction_T2_only"]] = 0
            singlePair["num_junction_shared"] = len(
                sorted_junction_df[sorted_junction_df["transcript_id"] == tx1_name]
            )
            singlePair["prop_junction_diff"] = 0
            singlePair["prop_junction_similar"] = 1
            singlePair[["junction_T1_only", "junction_T2_only"]] = ""
            singlePair[
                    ["junction_shared", "junction_T1_all", "junction_T2_all"]
                ] = "|".join(
                    sorted_junction_df[
                            sorted_junction_df["transcript_id"] == tx1_name
                        ]["coords"]
            )

    # Transcript pair does not share all junctions
    else:
        # Check if both transcripts are monoexon but do not overlap (not fsm)
        if len(sorted_junction_df) == 0:
            singlePair[
                ["num_junction_T1_only", "num_junction_T2_only", "num_junction_shared"]
            ] = 0
            singlePair["prop_junction_diff"] = 0
            singlePair["prop_junction_similar"] = 1
            singlePair[
                ["junction_T1_only",
                 "junction_T2_only",
                 "junction_shared",
                 "junction_T1_all",
                 "junction_T2_all"]
            ] = ""
        # Check if one of the transcripts is monoexon (only junctions from one transcript and not
        # the other)
        elif sorted_junction_df["transcript_id"].nunique() == 1:
            # Only T1 is monoexon
            if tx1_name in sorted_junction_df["transcript_id"].unique():
                singlePair[["num_junction_T1_only", "num_junction_shared"]] = 0
                singlePair["prop_junction_similar"] = 0
                singlePair["prop_junction_diff"] = 1
                singlePair[
                    ["junction_T1_only", "junction_T1_all", "junction_shared"]
                ] = ""
                singlePair["num_junction_T2_only"] = len(sorted_junction_df)
                singlePair[
                    ["junction_T2_only", "junction_T2_all"]
                ] = "|".join(sorted_junction_df["coords"])
            # Only T2 is monoexon
            else:
                singlePair[["num_junction_T2_only", "num_junction_shared"]] = 0
                singlePair["prop_junction_similar"] = 0
                singlePair["prop_junction_diff"] = 1
                singlePair[
                    ["junction_T2_only", "junction_T2_all", "junction_shared"]
                ] = ""
                singlePair["num_junction_T1_only"] = len(sorted_junction_df)
                singlePair[
                    ["junction_T1_only", "junction_T1_all"]
                ] = "|".join(sorted_junction_df["coords"])
        # Both transcripts multi-exon
        else:
            t1Junc = sorted_junction_df[
                    sorted_junction_df["transcript_id"] == tx1_name
                ]["coords"]
            t2Junc = sorted_junction_df[
                    sorted_junction_df["transcript_id"] == tx2_name
                ]["coords"]
            singlePair["junction_T1_all"] = "|".join(t1Junc)
            singlePair["junction_T2_all"] = "|".join(t2Junc)
            juncSharedSet = set(t1Junc).intersection(set(t2Junc))
            juncT1Set = set(t1Junc).difference(set(t2Junc))
            juncT2Set = set(t2Junc).difference(set(t1Junc))
            singlePair["num_junction_T1_only"] = len(juncT1Set)
            singlePair["num_junction_T2_only"] = len(juncT2Set)
            singlePair["num_junction_shared"] = len(juncSharedSet)
            singlePair["prop_junction_diff"] = (
                singlePair["num_junction_T1_only"] + singlePair["num_junction_T2_only"]
            ) / (
                singlePair["num_junction_T1_only"]
                + singlePair["num_junction_T2_only"]
                + singlePair["num_junction_shared"]
            )
            singlePair["prop_junction_similar"] = 1 - singlePair["prop_junction_diff"]
            singlePair["junction_T1_only"] = "|".join(
                sorted(juncT1Set, key=lambda x: t1Junc[t1Junc == x].index[0])
            )
            singlePair["junction_T2_only"] = "|".join(
                sorted(juncT2Set, key=lambda x: t2Junc[t2Junc == x].index[0])
            )
            singlePair["junction_shared"] = "|".join(
                sorted(juncSharedSet, key=lambda x: t2Junc[t2Junc == x].index[0])
            )
    return singlePair


def get_ER_distance(singlePair, ef_df, tx1_name, tx2_name, fsm=False):
    # Check if transcript pair shares all junctions (FSM, full-splice match)
    if fsm:
        # All junctions are shared, therefore all ER are shared
        singlePair[["num_ER_T1_only", "num_ER_T2_only"]] = 0
        singlePair["num_ER_shared"] = ef_df["er_id"].nunique()
        singlePair["prop_ER_diff"] = 0
        singlePair["prop_ER_similar"] = 1
        singlePair[["ER_T1_only", "ER_T2_only"]] = ""
        singlePair["ER_shared"] = "|".join(ef_df["region_id"].unique())
        ERSharedSet = ef_df["region_id"].drop_duplicates()
    else:
        # Check for shared and unique ER
        ERall = (
            ef_df.groupby("region_id")
            .agg({"transcript_id": (lambda x: max(x, key=len)), "flag_shared": "max"})
            .reset_index()
        )
        ERSharedSet = ERall[ERall["flag_shared"] == 1]["region_id"]
        ERT1Set = ERall[
            (ERall["flag_shared"] == 0) & (ERall["transcript_id"] == tx1_name)
        ]["region_id"]
        ERT2Set = ERall[
            (ERall["flag_shared"] == 0) & (ERall["transcript_id"] == tx2_name)
        ]["region_id"]
        singlePair["num_ER_T1_only"] = len(ERT1Set)
        singlePair["num_ER_T2_only"] = len(ERT2Set)
        singlePair["num_ER_shared"] = len(ERSharedSet)
        singlePair["prop_ER_diff"] = (
            singlePair["num_ER_T1_only"] + singlePair["num_ER_T2_only"]
        ) / (
            singlePair["num_ER_T1_only"]
            + singlePair["num_ER_T2_only"]
            + singlePair["num_ER_shared"]
        )
        singlePair["prop_ER_similar"] = 1 - singlePair["prop_ER_diff"]
        singlePair["ER_T1_only"] = "|".join(ERT1Set)
        singlePair["ER_T2_only"] = "|".join(ERT2Set)
        singlePair["ER_shared"] = "|".join(ERSharedSet)
    return singlePair, ERSharedSet


def get_EF_distance(singlePair, ef_df, tx1_name, tx2_name, ERSharedSet):
    num_shared = ef_df["flag_shared"].sum()
    # Check if all fragments are shared (transcripts are identical)
    if num_shared == len(ef_df):
        # All EF shared, transcripts are identical
        singlePair[["num_fragment_T1_only", "num_fragment_T2_only"]] = 0
        singlePair["num_fragment_shared"] = len(ef_df)
        singlePair["prop_fragment_diff"] = 0
        singlePair["prop_fragment_similar"] = 1
        singlePair[["fragment_T1_only", "fragment_T2_only"]] = ""
        singlePair["fragment_shared"] = "|".join(ef_df["fragment_id"])
        singlePair[
            ["num_fragment_singletons_T1_only", "num_fragment_singletons_T2_only"]
        ] = 0
        singlePair["num_fragment_singletons_shared"] = singlePair["num_fragment_shared"]

        # Nucleotide differences are 0, all are shared and are in shared ER
        singlePair[
            [
                "total_nt",
                "total_nt_in_shared_ER",
                "num_nt_shared",
                "num_nt_shared_in_shared_ER",
            ]
        ] = ef_df["fragment_length"].sum()
        singlePair[
            [
                "num_nt_T1_only",
                "num_nt_T1_only_in_shared_ER",
                "num_nt_T1_only_in_unique_ER",
                "num_nt_T2_only",
                "num_nt_T2_only_in_shared_ER",
                "num_nt_T2_only_in_unique_ER",
                "num_nt_diff",
            ]
        ] = 0
        singlePair[["prop_nt_diff", "prop_nt_diff_in_shared_ER"]] = 0
        singlePair[["prop_nt_similar", "prop_nt_similar_in_shared_ER"]] = 1

        # No IR present between identical transcripts
        singlePair[["num_IR_fragment_T1", "num_IR_fragment_T2"]] = 0
        singlePair[["IR_fragment_T1", "IR_fragment_T2"]] = ""
    else:
        # Not all EF shared
        fragSharedSet = ef_df[ef_df["flag_shared"] == 1]
        fragSharedSingSet = ef_df[
            (ef_df["flag_shared"] == 1) & (ef_df["flag_singleton"] == 1)
        ]
        fragT1Set = ef_df[
            (ef_df["flag_shared"] == 0) & (ef_df["transcript_id"] == tx1_name)
        ]
        fragT1SingSet = ef_df[
            (ef_df["flag_shared"] == 0)
            & (ef_df["transcript_id"] == tx1_name)
            & (ef_df["flag_singleton"] == 1)
        ]
        fragT2Set = ef_df[
            (ef_df["flag_shared"] == 0) & (ef_df["transcript_id"] == tx2_name)
        ]
        fragT2SingSet = ef_df[
            (ef_df["flag_shared"] == 0)
            & (ef_df["transcript_id"] == tx2_name)
            & (ef_df["flag_singleton"] == 1)
        ]
        singlePair["num_fragment_T1_only"] = len(fragT1Set)
        singlePair["num_fragment_T2_only"] = len(fragT2Set)
        singlePair["num_fragment_shared"] = len(fragSharedSet)
        singlePair["prop_fragment_diff"] = (
            singlePair["num_fragment_T1_only"] + singlePair["num_fragment_T2_only"]
        ) / (
            singlePair["num_fragment_T1_only"]
            + singlePair["num_fragment_T2_only"]
            + singlePair["num_fragment_shared"]
        )
        singlePair["prop_fragment_similar"] = 1 - singlePair["prop_fragment_diff"]
        singlePair["fragment_T1_only"] = "|".join(fragT1Set["fragment_id"])
        singlePair["fragment_T2_only"] = "|".join(fragT2Set["fragment_id"])
        singlePair["fragment_shared"] = "|".join(fragSharedSet["fragment_id"])
        singlePair["num_fragment_singletons_T1_only"] = len(fragT1SingSet)
        singlePair["num_fragment_singletons_T2_only"] = len(fragT2SingSet)
        singlePair["num_fragment_singletons_shared"] = len(fragSharedSingSet)

        # Count number of nt shared/different in all EF
        singlePair["num_nt_shared"] = fragSharedSet["fragment_length"].sum()
        singlePair["num_nt_T1_only"] = fragT1Set["fragment_length"].sum()
        singlePair["num_nt_T2_only"] = fragT2Set["fragment_length"].sum()
        singlePair["num_nt_diff"] = (
            singlePair["num_nt_T1_only"] + singlePair["num_nt_T2_only"]
        )
        singlePair["total_nt"] = (
            singlePair["num_nt_shared"]
            + singlePair["num_nt_T1_only"]
            + singlePair["num_nt_T2_only"]
        )
        singlePair["prop_nt_diff"] = (
            singlePair["num_nt_T1_only"] + singlePair["num_nt_T2_only"]
        ) / (singlePair["total_nt"])
        singlePair["prop_nt_similar"] = 1 - singlePair["prop_nt_diff"]

        # Count number of nt shared/different in EF only in shared ER
        singlePair["num_nt_T1_only_in_shared_ER"] = fragT1Set[
            fragT1Set["region_id"].isin(ERSharedSet)
        ]["fragment_length"].sum()
        singlePair["num_nt_T2_only_in_shared_ER"] = fragT2Set[
            fragT2Set["region_id"].isin(ERSharedSet)
        ]["fragment_length"].sum()
        singlePair["num_nt_shared_in_shared_ER"] = fragSharedSet[
            fragSharedSet["region_id"].isin(ERSharedSet)
        ]["fragment_length"].sum()
        singlePair["total_nt_in_shared_ER"] = (
            singlePair["num_nt_shared_in_shared_ER"]
            + singlePair["num_nt_T1_only_in_shared_ER"]
            + singlePair["num_nt_T2_only_in_shared_ER"]
        )
        if singlePair["total_nt_in_shared_ER"] != 0:
            singlePair["prop_nt_diff_in_shared_ER"] = (
                singlePair["num_nt_T1_only_in_shared_ER"]
                + singlePair["num_nt_T2_only_in_shared_ER"]
            ) / (singlePair["total_nt_in_shared_ER"])
            singlePair["prop_nt_similar_in_shared_ER"] = (
                1 - singlePair["prop_nt_diff_in_shared_ER"]
            )
        else:
            singlePair["prop_nt_diff_in_shared_ER"] = 0
            singlePair["prop_nt_similar_in_shared_ER"] = 0
        singlePair["num_nt_T1_only_in_unique_ER"] = fragT1Set[
            ~fragT1Set["region_id"].isin(ERSharedSet)
        ]["fragment_length"].sum()
        singlePair["num_nt_T2_only_in_unique_ER"] = fragT2Set[
            ~fragT2Set["region_id"].isin(ERSharedSet)
        ]["fragment_length"].sum()

        # Get IR distance values
        singlePair["num_IR_fragment_T1"] = len(
            fragT1Set[fragT1Set["ef_ir_flag"].map(int) == 1]
        )
        singlePair["num_IR_fragment_T2"] = len(
            fragT2Set[fragT2Set["ef_ir_flag"].map(int) == 1]
        )
        singlePair["IR_fragment_T1"] = "|".join(
            fragT1Set[fragT1Set["ef_ir_flag"].map(int) == 1]["fragment_id"]
        )
        singlePair["IR_fragment_T2"] = "|".join(
            fragT2Set[fragT2Set["ef_ir_flag"].map(int) == 1]["fragment_id"]
        )
    return singlePair


def set_AS_flags(singlePair, ef_df, tx1_name, tx2_name, has_junctions):
    # Flag alternative splicing events

    # Flag transcripts with no shared nucleotides (nonoverlapping)
    singlePair["flag_no_shared_nt"] = np.where(singlePair["prop_nt_diff"] == 1, 1, 0)

    # Flag full-splice matches (FSM, share all junctions)
    # NOTE: monoexon transcripts can be FSM
    singlePair["flag_FSM"] = np.where(singlePair["prop_junction_diff"] == 0, 1, 0)

    # Flag incomplete-splice matches
    # (ISM, one set of junctions is a complete consecutive subset of the other)
    # NOTE: ISM flags are 0 if transcripts are FSM or if at least one is monoexon
    if has_junctions:
        if (
            singlePair["flag_FSM"] == 1
            or singlePair["junction_T1_all"] == ""
            or singlePair["junction_T2_all"] == ""
        ):
            singlePair[
                ["flag_T1_ISM_of_T2", "flag_T2_ISM_of_T1"]
            ] = 0
        else:
            singlePair["flag_T1_ISM_of_T2"] = np.where(
                singlePair["junction_T1_all"] in singlePair["junction_T2_all"],
                1,
                0
            )
            singlePair["flag_T2_ISM_of_T1"] = np.where(
                singlePair["junction_T2_all"] in singlePair["junction_T1_all"],
                1,
                0
            )

    # If a pair contains at least one IR event - 5' and 3' variation calculated
    #   but flag_alt_exon and flag_alt_donor_acceptor set to 0
    # Intron retention
    singlePair["flag_IR"] = np.where(
        singlePair["num_IR_fragment_T1"] + singlePair["num_IR_fragment_T2"] > 0, 1, 0
    )

    # 5' and 3' end variation
    t1_min_start = ef_df[ef_df["transcript_id"].str.contains(tx1_name)][
        "ef_start"
    ].min()
    t2_min_start = ef_df[ef_df["transcript_id"].str.contains(tx2_name)][
        "ef_start"
    ].min()
    t1_max_end = ef_df[ef_df["transcript_id"].str.contains(tx1_name)]["ef_end"].max()
    t2_max_end = ef_df[ef_df["transcript_id"].str.contains(tx2_name)]["ef_end"].max()
    ef_strand = ef_df["ef_strand"].unique()[0]
    # 5' variation: where difference in start if + strand or difference in end if - strand
    singlePair["flag_5_variation"] = np.where(
        ((ef_strand == "+") & (t1_min_start != t2_min_start))
        | ((ef_strand == "-") & (t1_max_end != t2_max_end)),
        1,
        0,
    )
    # 3' end variation: where difference in end if + strand or difference in start if - strand
    singlePair["flag_3_variation"] = np.where(
        ((ef_strand == "-") & (t1_min_start != t2_min_start))
        | ((ef_strand == "+") & (t1_max_end != t2_max_end)),
        1,
        0,
    )

    # If a pair contains at least one IR event - 5' and 3' variation calculated
    #   but flag_alt_exon and flag_alt_donor_acceptor set to 0
    if singlePair["flag_IR"] == 1:
        singlePair["flag_alt_exon"] = 0
        singlePair["flag_alt_donor_acceptor"] = 0
    else:
        # Alternate exons are when not all exon regions are shared
        singlePair["flag_alt_exon"] = np.where(singlePair["prop_ER_diff"] > 0, 1, 0)

        # Alternate donor/acceptors in shared ER (where not all junctions are shared)
        singlePair["flag_alt_donor_acceptor"] = np.where(
            (singlePair["prop_nt_diff_in_shared_ER"] > 0)
            & (singlePair["flag_FSM"] == 0),
            1,
            0,
        )

    # If transcripts are nonoverlapping then 5'/3' variation and alt exon flags set to 0
    if singlePair["flag_no_shared_nt"] == 1:
        singlePair["flag_alt_exon"] = 0
        singlePair["flag_5_variation"] = 0
        singlePair["flag_3_variation"] = 0
        singlePair["flag_FSM"] = 0
        singlePair["flag_T1_ISM_of_T2"] = 0
        singlePair["flag_T2_ISM_of_T1"] = 0

    return singlePair

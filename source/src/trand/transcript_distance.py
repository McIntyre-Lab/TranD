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
    "num_jxn_only_T1",
    "num_jxn_only_T2",
    "num_jxn_ovlp",
    "prop_jxn_noOvlp",
    "prop_jxn_ovlp",
    "jxn_string_T1",
    "jxn_string_T2",
    "jxn_only_T1",
    "jxn_only_T2",
    "jxn_same",
    "num_ER_only_T1",
    "num_ER_only_T2",
    "num_ER_ovlp",
    "prop_ER_noOvlp",
    "prop_ER_ovlp",
    "ER_only_T1",
    "ER_only_T2",
    "ER_ovlp",
    "num_EF_only_T1",
    "num_EF_only_T2",
    "num_EF_ovlp",
    "prop_EF_noOvlp",
    "prop_EF_ovlp",
    "EF_only_T1",
    "EF_only_T2",
    "EF_ovlp",
    "num_exon_only_T1",
    "num_exon_only_T2",
    "num_exon_ovlp",
    "num_IR_EF_T1",
    "num_IR_EF_T2",
    "IR_EF_T1",
    "IR_EF_T2",
    "num_nt_ovlp",
    "num_nt_only_T1",
    "num_nt_only_T2",
    "num_nt_noOvlp",
    "total_nt",
    "prop_nt_noOvlp",
    "prop_nt_ovlp",
    "num_nt_only_T1_in_ovlpER",
    "num_nt_only_T2_in_ovlpER",
    "num_nt_ovlp_in_ovlpER",
    "total_nt_in_ovlpER",
    "prop_nt_noOvlp_in_ovlpER",
    "prop_nt_ovlp_in_ovlpER",
    "num_nt_only_T1_in_uniqER",
    "num_nt_only_T2_in_uniqER",
    "flag_FSM",
    "flag_T1_ISM_of_T2",
    "flag_T2_ISM_of_T1",
    "flag_IR",
    "flag_5_var",
    "flag_3_var",
    "flag_alt_DA",
    "flag_alt_exon",
    "flag_no_ovlp_nt",
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

    # Set exon EF and exon region ID values based on genomic coordinates
    ef_df["EF_id"] = (
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

    # Flag shared EFs
    ef_df["flag_ovlp"] = np.where(
        (ef_df["transcript_id"] == tx1_name + "|" + tx2_name)
        | (ef_df["transcript_id"] == tx2_name + "|" + tx1_name),
        1,
        0,
    )

    # Flag singleton EFs (where the EF is a full exon region)
    ef_df["num_EF_in_ER"] = ef_df.groupby("er_id")["ef_id"].transform("count")
    ef_df["flag_singleton"] = np.where(ef_df["num_EF_in_ER"] == 1, 1, 0)

    # Add EF lengths (end - start)
    ef_df["EF_length"] = ef_df["ef_end"].map(int) - ef_df["ef_start"].map(int)

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

    # Gt distance measures for exon regions (ER) and exon EFs (EF)
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
                ["num_jxn_only_T1", "num_jxn_only_T2", "num_jxn_ovlp"]
            ] = 0
            singlePair["prop_jxn_noOvlp"] = 0
            singlePair["prop_jxn_ovlp"] = 1
            singlePair[
                ["jxn_only_T1",
                 "jxn_only_T2",
                 "jxn_same",
                 "jxn_string_T1",
                 "jxn_string_T2"]
            ] = ""
        # FSM transcripts are multiexon
        else:
            singlePair[["num_jxn_only_T1", "num_jxn_only_T2"]] = 0
            singlePair["num_jxn_ovlp"] = len(
                sorted_junction_df[sorted_junction_df["transcript_id"] == tx1_name]
            )
            singlePair["prop_jxn_noOvlp"] = 0
            singlePair["prop_jxn_ovlp"] = 1
            singlePair[["jxn_only_T1", "jxn_only_T2"]] = ""
            singlePair[
                    ["jxn_same", "jxn_string_T1", "jxn_string_T2"]
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
                ["num_jxn_only_T1", "num_jxn_only_T2", "num_jxn_ovlp"]
            ] = 0
            singlePair["prop_jxn_noOvlp"] = 0
            singlePair["prop_jxn_ovlp"] = 1
            singlePair[
                ["jxn_only_T1",
                 "jxn_only_T2",
                 "jxn_same",
                 "jxn_string_T1",
                 "jxn_string_T2"]
            ] = ""
        # Check if one of the transcripts is monoexon (only junctions from one transcript and not
        # the other)
        elif sorted_junction_df["transcript_id"].nunique() == 1:
            # Only T1 is monoexon
            if tx2_name in sorted_junction_df["transcript_id"].unique():
                singlePair[["num_jxn_only_T1", "num_jxn_ovlp"]] = 0
                singlePair["prop_jxn_ovlp"] = 0
                singlePair["prop_jxn_noOvlp"] = 1
                singlePair[
                    ["jxn_only_T1", "jxn_string_T1", "jxn_same"]
                ] = ""
                singlePair["num_jxn_only_T2"] = len(sorted_junction_df)
                singlePair[
                    ["jxn_only_T2", "jxn_string_T2"]
                ] = "|".join(sorted_junction_df["coords"])
            # Only T2 is monoexon
            else:
                singlePair[["num_jxn_only_T2", "num_jxn_ovlp"]] = 0
                singlePair["prop_jxn_ovlp"] = 0
                singlePair["prop_jxn_noOvlp"] = 1
                singlePair[
                    ["jxn_only_T2", "jxn_string_T2", "jxn_same"]
                ] = ""
                singlePair["num_jxn_only_T1"] = len(sorted_junction_df)
                singlePair[
                    ["jxn_only_T1", "jxn_string_T1"]
                ] = "|".join(sorted_junction_df["coords"])
        # Both transcripts multi-exon
        else:
            t1Junc = sorted_junction_df[
                    sorted_junction_df["transcript_id"] == tx1_name
                ]["coords"]
            t2Junc = sorted_junction_df[
                    sorted_junction_df["transcript_id"] == tx2_name
                ]["coords"]
            singlePair["jxn_string_T1"] = "|".join(t1Junc)
            singlePair["jxn_string_T2"] = "|".join(t2Junc)
            juncSharedSet = set(t1Junc).intersection(set(t2Junc))
            juncT1Set = set(t1Junc).difference(set(t2Junc))
            juncT2Set = set(t2Junc).difference(set(t1Junc))
            singlePair["num_jxn_only_T1"] = len(juncT1Set)
            singlePair["num_jxn_only_T2"] = len(juncT2Set)
            singlePair["num_jxn_ovlp"] = len(juncSharedSet)
            singlePair["prop_jxn_noOvlp"] = (
                singlePair["num_jxn_only_T1"] + singlePair["num_jxn_only_T2"]
            ) / (
                singlePair["num_jxn_only_T1"]
                + singlePair["num_jxn_only_T2"]
                + singlePair["num_jxn_ovlp"]
            )
            singlePair["prop_jxn_ovlp"] = 1 - singlePair["prop_jxn_noOvlp"]
            singlePair["jxn_only_T1"] = "|".join(
                sorted(juncT1Set, key=lambda x: t1Junc[t1Junc == x].index[0])
            )
            singlePair["jxn_only_T2"] = "|".join(
                sorted(juncT2Set, key=lambda x: t2Junc[t2Junc == x].index[0])
            )
            singlePair["jxn_same"] = "|".join(
                sorted(juncSharedSet, key=lambda x: t2Junc[t2Junc == x].index[0])
            )
    return singlePair


def get_ER_distance(singlePair, ef_df, tx1_name, tx2_name, fsm=False):
    # Check if transcript pair shares all junctions (FSM, full-splice match)
    if fsm:
        # All junctions are shared, therefore all ER are shared
        singlePair[["num_ER_only_T1", "num_ER_only_T2"]] = 0
        singlePair["num_ER_ovlp"] = ef_df["er_id"].nunique()
        singlePair["prop_ER_noOvlp"] = 0
        singlePair["prop_ER_ovlp"] = 1
        singlePair[["ER_only_T1", "ER_only_T2"]] = ""
        singlePair["ER_ovlp"] = "|".join(ef_df["region_id"].unique())
        ERSharedSet = ef_df["region_id"].drop_duplicates()
    else:
        # Check for shared and unique ER
        ERall = (
            ef_df.groupby("region_id")
            .agg({"transcript_id": (lambda x: max(x, key=len)), "flag_ovlp": "max"})
            .reset_index()
        )
        ERSharedSet = ERall[ERall["flag_ovlp"] == 1]["region_id"]
        ERT1Set = ERall[
            (ERall["flag_ovlp"] == 0) & (ERall["transcript_id"] == tx1_name)
        ]["region_id"]
        ERT2Set = ERall[
            (ERall["flag_ovlp"] == 0) & (ERall["transcript_id"] == tx2_name)
        ]["region_id"]
        singlePair["num_ER_only_T1"] = len(ERT1Set)
        singlePair["num_ER_only_T2"] = len(ERT2Set)
        singlePair["num_ER_ovlp"] = len(ERSharedSet)
        singlePair["prop_ER_noOvlp"] = (
            singlePair["num_ER_only_T1"] + singlePair["num_ER_only_T2"]
        ) / (
            singlePair["num_ER_only_T1"]
            + singlePair["num_ER_only_T2"]
            + singlePair["num_ER_ovlp"]
        )
        singlePair["prop_ER_ovlp"] = 1 - singlePair["prop_ER_noOvlp"]
        singlePair["ER_only_T1"] = "|".join(ERT1Set)
        singlePair["ER_only_T2"] = "|".join(ERT2Set)
        singlePair["ER_ovlp"] = "|".join(ERSharedSet)
    return singlePair, ERSharedSet


def get_EF_distance(singlePair, ef_df, tx1_name, tx2_name, ERSharedSet):
    num_ovlp = ef_df["flag_ovlp"].sum()
    # Check if all EFs are shared (transcripts are identical)
    if num_ovlp == len(ef_df):
        # All EF shared, transcripts are identical
        singlePair[["num_EF_only_T1", "num_EF_only_T2"]] = 0
        singlePair["num_EF_ovlp"] = len(ef_df)
        singlePair["prop_EF_noOvlp"] = 0
        singlePair["prop_EF_ovlp"] = 1
        singlePair[["EF_only_T1", "EF_only_T2"]] = ""
        singlePair["EF_ovlp"] = "|".join(ef_df["EF_id"])
        singlePair[
            ["num_exon_only_T1", "num_exon_only_T2"]
        ] = 0
        singlePair["num_exon_ovlp"] = singlePair["num_EF_ovlp"]

        # Nucleotide differences are 0, all are shared and are in shared ER
        singlePair[
            [
                "total_nt",
                "total_nt_in_ovlpER",
                "num_nt_ovlp",
                "num_nt_ovlp_in_ovlpER",
            ]
        ] = ef_df["EF_length"].sum()
        singlePair[
            [
                "num_nt_only_T1",
                "num_nt_only_T1_in_ovlpER",
                "num_nt_only_T1_in_uniqER",
                "num_nt_only_T2",
                "num_nt_only_T2_in_ovlpER",
                "num_nt_only_T2_in_uniqER",
                "num_nt_noOvlp",
            ]
        ] = 0
        singlePair[["prop_nt_noOvlp", "prop_nt_noOvlp_in_ovlpER"]] = 0
        singlePair[["prop_nt_ovlp", "prop_nt_ovlp_in_ovlpER"]] = 1

        # No IR present between identical transcripts
        singlePair[["num_IR_EF_T1", "num_IR_EF_T2"]] = 0
        singlePair[["IR_EF_T1", "IR_EF_T2"]] = ""
    else:
        # Not all EF shared
        fragSharedSet = ef_df[ef_df["flag_ovlp"] == 1]
        fragSharedSingSet = ef_df[
            (ef_df["flag_ovlp"] == 1) & (ef_df["flag_singleton"] == 1)
        ]
        fragT1Set = ef_df[
            (ef_df["flag_ovlp"] == 0) & (ef_df["transcript_id"] == tx1_name)
        ]
        fragT1SingSet = ef_df[
            (ef_df["flag_ovlp"] == 0)
            & (ef_df["transcript_id"] == tx1_name)
            & (ef_df["flag_singleton"] == 1)
        ]
        fragT2Set = ef_df[
            (ef_df["flag_ovlp"] == 0) & (ef_df["transcript_id"] == tx2_name)
        ]
        fragT2SingSet = ef_df[
            (ef_df["flag_ovlp"] == 0)
            & (ef_df["transcript_id"] == tx2_name)
            & (ef_df["flag_singleton"] == 1)
        ]
        singlePair["num_EF_only_T1"] = len(fragT1Set)
        singlePair["num_EF_only_T2"] = len(fragT2Set)
        singlePair["num_EF_ovlp"] = len(fragSharedSet)
        singlePair["prop_EF_noOvlp"] = (
            singlePair["num_EF_only_T1"] + singlePair["num_EF_only_T2"]
        ) / (
            singlePair["num_EF_only_T1"]
            + singlePair["num_EF_only_T2"]
            + singlePair["num_EF_ovlp"]
        )
        singlePair["prop_EF_ovlp"] = 1 - singlePair["prop_EF_noOvlp"]
        singlePair["EF_only_T1"] = "|".join(fragT1Set["EF_id"])
        singlePair["EF_only_T2"] = "|".join(fragT2Set["EF_id"])
        singlePair["EF_ovlp"] = "|".join(fragSharedSet["EF_id"])
        singlePair["num_exon_only_T1"] = len(fragT1SingSet)
        singlePair["num_exon_only_T2"] = len(fragT2SingSet)
        singlePair["num_exon_ovlp"] = len(fragSharedSingSet)

        # Count number of nt shared/different in all EF
        singlePair["num_nt_ovlp"] = fragSharedSet["EF_length"].sum()
        singlePair["num_nt_only_T1"] = fragT1Set["EF_length"].sum()
        singlePair["num_nt_only_T2"] = fragT2Set["EF_length"].sum()
        singlePair["num_nt_noOvlp"] = (
            singlePair["num_nt_only_T1"] + singlePair["num_nt_only_T2"]
        )
        singlePair["total_nt"] = (
            singlePair["num_nt_ovlp"]
            + singlePair["num_nt_only_T1"]
            + singlePair["num_nt_only_T2"]
        )
        singlePair["prop_nt_noOvlp"] = (
            singlePair["num_nt_only_T1"] + singlePair["num_nt_only_T2"]
        ) / (singlePair["total_nt"])
        singlePair["prop_nt_ovlp"] = 1 - singlePair["prop_nt_noOvlp"]

        # Count number of nt shared/different in EF only in shared ER
        singlePair["num_nt_only_T1_in_ovlpER"] = fragT1Set[
            fragT1Set["region_id"].isin(ERSharedSet)
        ]["EF_length"].sum()
        singlePair["num_nt_only_T2_in_ovlpER"] = fragT2Set[
            fragT2Set["region_id"].isin(ERSharedSet)
        ]["EF_length"].sum()
        singlePair["num_nt_ovlp_in_ovlpER"] = fragSharedSet[
            fragSharedSet["region_id"].isin(ERSharedSet)
        ]["EF_length"].sum()
        singlePair["total_nt_in_ovlpER"] = (
            singlePair["num_nt_ovlp_in_ovlpER"]
            + singlePair["num_nt_only_T1_in_ovlpER"]
            + singlePair["num_nt_only_T2_in_ovlpER"]
        )
        if singlePair["total_nt_in_ovlpER"] != 0:
            singlePair["prop_nt_noOvlp_in_ovlpER"] = (
                singlePair["num_nt_only_T1_in_ovlpER"]
                + singlePair["num_nt_only_T2_in_ovlpER"]
            ) / (singlePair["total_nt_in_ovlpER"])
            singlePair["prop_nt_ovlp_in_ovlpER"] = (
                1 - singlePair["prop_nt_noOvlp_in_ovlpER"]
            )
        else:
            singlePair["prop_nt_noOvlp_in_ovlpER"] = 0
            singlePair["prop_nt_ovlp_in_ovlpER"] = 0
        singlePair["num_nt_only_T1_in_uniqER"] = fragT1Set[
            ~fragT1Set["region_id"].isin(ERSharedSet)
        ]["EF_length"].sum()
        singlePair["num_nt_only_T2_in_uniqER"] = fragT2Set[
            ~fragT2Set["region_id"].isin(ERSharedSet)
        ]["EF_length"].sum()

        # Get IR distance values
        singlePair["num_IR_EF_T1"] = len(
            fragT1Set[fragT1Set["ef_ir_flag"].map(int) == 1]
        )
        singlePair["num_IR_EF_T2"] = len(
            fragT2Set[fragT2Set["ef_ir_flag"].map(int) == 1]
        )
        singlePair["IR_EF_T1"] = "|".join(
            fragT1Set[fragT1Set["ef_ir_flag"].map(int) == 1]["EF_id"]
        )
        singlePair["IR_EF_T2"] = "|".join(
            fragT2Set[fragT2Set["ef_ir_flag"].map(int) == 1]["EF_id"]
        )
    return singlePair


def set_AS_flags(singlePair, ef_df, tx1_name, tx2_name, has_junctions):
    # Flag alternative splicing events

    # Flag transcripts with no shared nucleotides (nonoverlapping)
    singlePair["flag_no_ovlp_nt"] = np.where(singlePair["prop_nt_noOvlp"] == 1, 1, 0)

    # Flag full-splice matches (FSM, share all junctions)
    # NOTE: monoexon transcripts can be FSM
    singlePair["flag_FSM"] = np.where(singlePair["prop_jxn_noOvlp"] == 0, 1, 0)

    # Flag incomplete-splice matches
    # (ISM, one set of junctions is a complete consecutive subset of the other)
    # NOTE: ISM flags are 0 if transcripts are FSM or if at least one is monoexon
    if has_junctions:
        if (
            singlePair["flag_FSM"] == 1
            or singlePair["jxn_string_T1"] == ""
            or singlePair["jxn_string_T2"] == ""
        ):
            singlePair[
                ["flag_T1_ISM_of_T2", "flag_T2_ISM_of_T1"]
            ] = 0
        else:
            singlePair["flag_T1_ISM_of_T2"] = np.where(
                singlePair["jxn_string_T1"] in singlePair["jxn_string_T2"],
                1,
                0
            )
            singlePair["flag_T2_ISM_of_T1"] = np.where(
                singlePair["jxn_string_T2"] in singlePair["jxn_string_T1"],
                1,
                0
            )

    # If a pair contains at least one IR event - 5' and 3' variation calculated
    #   but flag_alt_exon and flag_alt_DA set to 0
    # Intron retention
    singlePair["flag_IR"] = np.where(
        singlePair["num_IR_EF_T1"] + singlePair["num_IR_EF_T2"] > 0, 1, 0
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
    singlePair["flag_5_var"] = np.where(
        ((ef_strand == "+") & (t1_min_start != t2_min_start))
        | ((ef_strand == "-") & (t1_max_end != t2_max_end)),
        1,
        0,
    )
    # 3' end variation: where difference in end if + strand or difference in start if - strand
    singlePair["flag_3_var"] = np.where(
        ((ef_strand == "-") & (t1_min_start != t2_min_start))
        | ((ef_strand == "+") & (t1_max_end != t2_max_end)),
        1,
        0,
    )

    # If a pair contains at least one IR event - 5' and 3' variation calculated
    #   but flag_alt_exon and flag_alt_DA set to 0
    if singlePair["flag_IR"] == 1:
        singlePair["flag_alt_exon"] = 0
        singlePair["flag_alt_DA"] = 0
    else:
        # Alternate exons are when not all exon regions are shared
        singlePair["flag_alt_exon"] = np.where(singlePair["prop_ER_noOvlp"] > 0, 1, 0)

        # Alternate donor/acceptors in shared ER (where not all junctions are shared)
        singlePair["flag_alt_DA"] = np.where(
            (singlePair["prop_nt_noOvlp_in_ovlpER"] > 0)
            & (singlePair["flag_FSM"] == 0),
            1,
            0,
        )

    # If transcripts are nonoverlapping then 5'/3' variation and alt exon flags set to 0
    if singlePair["flag_no_ovlp_nt"] == 1:
        singlePair["flag_alt_exon"] = 0
        singlePair["flag_5_var"] = 0
        singlePair["flag_3_var"] = 0
        singlePair["flag_FSM"] = 0
        singlePair["flag_T1_ISM_of_T2"] = 0
        singlePair["flag_T2_ISM_of_T1"] = 0

    return singlePair

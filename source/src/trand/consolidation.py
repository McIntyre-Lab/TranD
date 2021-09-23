#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:38:56 2021

@author: Adalena Nanni <adalena.nanni@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu
"""

import pandas as pd
import numpy as np
from loguru import logger
from .io import open_output_files
from .io import write_gtf
from .io import write_output
from .bedtools import prep_bed_for_ea

# CONFIGURATION
consol_key_cols = ['gene_id', 'transcript_id', 'consolidation_transcript_id']

# Duplicate to prevent a circular import
jct_df_cols = ['gene_id', 'transcript_id', 'coords']


def create_junction_catalog_str(gene, tx, tx_data):
    """Create a junction catalog for a transcript"""
    junctions = []
    id = 0
    tx_data_df = pd.DataFrame(tx_data, columns=['chrom', 'start', 'end', 'name',
                              'score', 'strand'])
    tx_data_df['start'] = tx_data_df['start'].astype(int)
    tx_data_df['end'] = tx_data_df['end'].astype(int)
    tx_data_df = tx_data_df.sort_values('start')
    for index, e in tx_data_df.iterrows():
        id += 1
        if id == 1:
            left_end = int(e.end) - 10
            continue
        right_start = int(e.start + 10)
        jct_coords = f"{e.chrom}:{left_end}:{right_start}:{e.strand}"
        junctions.append([gene, tx, jct_coords])
        left_end = int(e.end) - 10
    return junctions


def consolidate_junctions(
    bed_gene_data, pre_consol_jct_df, outdir, skip_interm, prefix="tr"
):
    gene_id = bed_gene_data["gene_id"]
    logger.debug("Consolidation of {}", gene_id)

    # Set output data dataframe
    consol_gene_cols = ["seqname", "start", "end", "strand", "gene_id", "transcript_id"]
    consol_gene = pd.DataFrame(columns=consol_gene_cols)

    # Check if single transcript, if yes then it cannot be consolidated
    if len(bed_gene_data["transcript_list"]) == 1:
        logger.debug("Gene {} has a single transcript.", gene_id)
        xcrpt_all = pd.DataFrame(
            bed_gene_data["transcript_list"], columns=["transcript_id"]
        )
        xcrpt_all["gene_id"] = gene_id
        xcrpt_all["chr"] = bed_gene_data[xcrpt_all["transcript_id"][0]][0][0]
        xcrpt_all["strand"] = bed_gene_data[xcrpt_all["transcript_id"][0]][0][5]
        xcrpt_all["consolidation_transcript_id"] = prefix + "_" + gene_id + "_1"
        key_gene = xcrpt_all[
            ["gene_id", "transcript_id", "consolidation_transcript_id"]
        ]
        # Format output data
        for consol_transcript, transcript in (
            key_gene.groupby("consolidation_transcript_id")["transcript_id"]
            .first()
            .iteritems()
        ):
            for row in bed_gene_data[transcript]:
                consol_gene = pd.concat(
                    [
                        consol_gene,
                        pd.DataFrame(
                            [
                                [
                                    row[0],
                                    int(row[1]) + 1,
                                    int(row[2]),
                                    row[5],
                                    gene_id,
                                    consol_transcript,
                                ]
                            ],
                            columns=consol_gene.columns,
                        ),
                    ],
                    ignore_index=True,
                )
                consol_gene = consol_gene.sort_values(
                    by=["transcript_id", "start"]
                ).reset_index(drop=True)

    # Check if more than one transcript present and all are monoexon (no junctions present)
    elif len(bed_gene_data["transcript_list"]) > 1 and len(pre_consol_jct_df) == 0:
        logger.debug(
            "Gene {} has all monoexon transcripts (more than 1 transcript present).",
            gene_id,
        )
        # Get transcript list and start/end of transcripts
        xcrpt_all = pd.DataFrame(
            bed_gene_data["transcript_list"], columns=["transcript_id"]
        )
        xcrpt_all["gene_id"] = gene_id
        xcrpt_all["chr"] = bed_gene_data[xcrpt_all["transcript_id"][0]][0][0]
        xcrpt_all["strand"] = bed_gene_data[xcrpt_all["transcript_id"][0]][0][5]
        xcrpt_all["start"] = xcrpt_all["transcript_id"].map(
            lambda xcrpt: min([int(x[1]) + 1 for x in bed_gene_data[xcrpt]])
        )
        xcrpt_all["end"] = xcrpt_all["transcript_id"].map(
            lambda xcrpt: max([int(x[2]) for x in bed_gene_data[xcrpt]])
        )
        # Loop over coordinates to make bins of overlapping monoexon transcripts
        xcrpt_all = xcrpt_all.sort_values(by=["start", "end"]).reset_index(drop=True)
        bin_num = 1
        bin_start = xcrpt_all["start"][0]
        bin_end = xcrpt_all["end"][0]
        xcrpt_all["monoexon_group"] = 0
        for index, row in xcrpt_all.iterrows():
            if row["start"] >= bin_start and row["start"] <= bin_end:
                xcrpt_all.loc[index, "monoexon_group"] = bin_num
                if row["end"] >= bin_end:
                    bin_end = row["end"]
            else:
                bin_num += 1
                bin_start = row["start"]
                bin_end = row["end"]
                xcrpt_all.loc[index, "monoexon_group"] = bin_num
        # Group by assigned bins to get min start and max end
        xcrpt_all["min_group_start"] = xcrpt_all.groupby("monoexon_group")[
            "start"
        ].transform("min")
        xcrpt_all["max_group_end"] = xcrpt_all.groupby("monoexon_group")[
            "end"
        ].transform("max")
        xcrpt_all["consolidation_transcript_id"] = (
            prefix
            + "_"
            + gene_id
            + "_"
            + xcrpt_all["monoexon_group"].astype(int).map(str)
        )
        key_gene = xcrpt_all[
            ["gene_id", "transcript_id", "consolidation_transcript_id"]
        ]
        # Format output data
        consol_xcrpt = (
            xcrpt_all.groupby("monoexon_group")[
                ["gene_id", "chr", "strand", "min_group_start", "max_group_end"]
            ]
            .first()
            .reset_index()
        )
        consol_xcrpt["consolidation_transcript_id"] = (
            prefix
            + "_"
            + gene_id
            + "_"
            + consol_xcrpt["monoexon_group"].astype(int).map(str)
        )
        consol_gene = consol_xcrpt[
            [
                "chr",
                "min_group_start",
                "max_group_end",
                "strand",
                "gene_id",
                "consolidation_transcript_id",
            ]
        ]
        consol_gene.columns = consol_gene_cols
    # Gene contains at least one multiexon transcript and more than one transcript total
    else:
        logger.debug(
            "Gene {} contains at least one multiexon transcript (more than 1 transcript present).",
            gene_id,
        )
        # Get junction string for each transcript in the gene
        pre_consol_jct_df = pre_consol_jct_df.sort_values(
            ["transcript_id", "coords"]
        ).reset_index(drop=True)
        collapse_pre_jct = (
            pre_consol_jct_df.groupby("transcript_id")
            .apply(func=lambda x: "|".join(x["coords"]))
            .reset_index()
            .rename(columns={0: "junctionID_order"})
            .sort_values(["junctionID_order"])
            .reset_index(drop=True)
        )

        # Merge transcript-level junction (multiexon transcripts only) and
        #     transcript_id values (all transcripts including multi-exon) variables to get
        #     mono-exon transcripts
        xcrpt_all = pd.DataFrame(
            bed_gene_data["transcript_list"], columns=["transcript_id"]
        )
        xcrpt_all["gene_id"] = gene_id
        xcrpt_all["chr"] = bed_gene_data[xcrpt_all["transcript_id"][0]][0][0]
        xcrpt_all["strand"] = bed_gene_data[xcrpt_all["transcript_id"][0]][0][5]
        xcrpt_jct_w_mono = pd.merge(
            xcrpt_all, collapse_pre_jct, how="outer", on="transcript_id", validate="1:1"
        )
        xcrpt_jct_w_mono["start"] = xcrpt_jct_w_mono["transcript_id"].map(
            lambda xcrpt: min([int(x[1]) + 1 for x in bed_gene_data[xcrpt]])
        )
        xcrpt_jct_w_mono["end"] = xcrpt_jct_w_mono["transcript_id"].map(
            lambda xcrpt: max([int(x[2]) for x in bed_gene_data[xcrpt]])
        )

        # Flag monoexon transcripts (no junctionID_order present)
        xcrpt_jct_w_mono["flag_monoexon_transcript"] = np.where(
            xcrpt_jct_w_mono["junctionID_order"].isna(), 1, 0
        )

        # If present, get monoexon junctionID_order values
        #   using group number for overlapping monoexon transcripts
        if xcrpt_jct_w_mono["flag_monoexon_transcript"].sum() > 0:
            # Loop over coordinates to make bins of overlapping monoexon transcripts
            xcrpt_jct_w_mono = xcrpt_jct_w_mono.sort_values(
                by=["start", "end"]
            ).reset_index(drop=True)
            bin_num = 1
            bin_start = xcrpt_jct_w_mono[
                xcrpt_jct_w_mono["flag_monoexon_transcript"] == 1
            ]["start"].values[0]
            bin_end = xcrpt_jct_w_mono[
                xcrpt_jct_w_mono["flag_monoexon_transcript"] == 1
            ]["end"].values[0]
            xcrpt_jct_w_mono["monoexon_group"] = 0
            for index, row in xcrpt_jct_w_mono[
                xcrpt_jct_w_mono["flag_monoexon_transcript"] == 1
            ].iterrows():
                if row["start"] >= bin_start and row["start"] <= bin_end:
                    xcrpt_jct_w_mono.loc[index, "monoexon_group"] = bin_num
                    if row["end"] >= bin_end:
                        bin_end = row["end"]
                else:
                    bin_num += 1
                    bin_start = row["start"]
                    bin_end = row["end"]
                    xcrpt_jct_w_mono.loc[index, "monoexon_group"] = bin_num
            logger.debug(
                "Gene {} has {} monoexon transcript(s) that occupy {} nonoverlapping "
                "genomic region(s).",
                gene_id,
                xcrpt_jct_w_mono["flag_monoexon_transcript"].sum(),
                bin_num,
            )
            # Set junctionID_order to monoexon group number
            xcrpt_jct_w_mono = xcrpt_jct_w_mono.fillna(-1)
            xcrpt_jct_w_mono.loc[
                xcrpt_jct_w_mono["junctionID_order"] == -1, "junctionID_order"
            ] = xcrpt_jct_w_mono["monoexon_group"]

        # Make groups of transcripts with same junctionID_order
        # Group names will be piped list of transcript IDs that share the same junctions
        xcrpt_jct_w_mono["transcriptID_cat"] = xcrpt_jct_w_mono.groupby(
            ["junctionID_order"]
        )["transcript_id"].transform(func=lambda x: "|".join(x))

        # Get minimum start and maximum end for longest 5' and 3' ends within each group
        #     with the same junctionID_order
        xcrpt_jct_w_mono["transcript_length"] = xcrpt_jct_w_mono["end"].astype(
            int
        ) - xcrpt_jct_w_mono["start"].astype(int)
        xcrpt_jct_w_mono["min_group_start"] = xcrpt_jct_w_mono.groupby(
            "junctionID_order"
        )["start"].transform("min")
        xcrpt_jct_w_mono["max_group_end"] = xcrpt_jct_w_mono.groupby(
            "junctionID_order"
        )["end"].transform("max")

        # Select the longest transcript for each junctionID_order
        longest_df = xcrpt_jct_w_mono.loc[
            xcrpt_jct_w_mono.groupby("junctionID_order")["transcript_length"].idxmax()
        ]

        # Generate new transcript_id for each transcript:
        #     [prefix]_[gene_id]_# where #={1,...,n} for all n transcripts in the gene
        # First sort transcript_length
        longest_df["consol_transcript_length"] = (
            longest_df["max_group_end"] - longest_df["min_group_start"]
        )
        longest_df = longest_df.sort_values(
            by="consol_transcript_length", ascending=False
        ).reset_index(drop=True)
        longest_df["transcript_rank_in_gene"] = longest_df[
            "consol_transcript_length"
        ].rank(method="first")
        longest_df["consolidation_transcript_id"] = (
            prefix
            + "_"
            + gene_id
            + "_"
            + longest_df["transcript_rank_in_gene"].astype(int).map(str)
        )

        # Flag transcripts where min start does not match the longest representative start

        longest_df["flag_not_min_start"] = np.where(
            longest_df["start"] != longest_df["min_group_start"], 1, 0
        )
        # Flag transcripts where the max end does not match the end of the longest representative
        longest_df["flag_not_max_end"] = np.where(
            longest_df["end"] != longest_df["max_group_end"], 1, 0
        )

        # Set consolidated start and end
        longest_df["consol_start"] = np.where(
            longest_df["flag_not_min_start"] == 1,
            longest_df["min_group_start"],
            longest_df["start"],
        )
        longest_df["consol_end"] = np.where(
            longest_df["flag_not_max_end"] == 1,
            longest_df["max_group_end"],
            longest_df["end"],
        )

        # Make key for gene_id 2 transcript_id 2 consolidation_transcript_id
        key_gene = pd.merge(
            xcrpt_jct_w_mono,
            longest_df[["consolidation_transcript_id", "junctionID_order"]],
            how="left",
            on="junctionID_order",
            validate="m:1",
        )
        key_gene["gene_id"] = gene_id
        key_gene = key_gene[["gene_id", "transcript_id", "consolidation_transcript_id"]]

        # Format output data by setting first and last exon coordinates based on longest first and
        # last exons of the groups Internal exons are the same for all transcripts in the group so
        # coordinates are kept for those.
        for consol_transcript, transcript in (
            key_gene.groupby("consolidation_transcript_id")["transcript_id"]
            .first()
            .iteritems()
        ):
            num = 0
            # Make start and end values integers for proper sorting
            for val in range(0, len(bed_gene_data[transcript])):
                tupleList = list(bed_gene_data[transcript][val])
                tupleList[1] = int(tupleList[1])
                tupleList[2] = int(tupleList[2])
                bed_gene_data[transcript][val] = tuple(tupleList)
            for row in sorted(bed_gene_data[transcript]):
                num = num + 1
                if num == 1:
                    consol_gene = pd.concat(
                        [
                            consol_gene,
                            pd.DataFrame(
                                [
                                    [
                                        row[0],
                                        longest_df[
                                            longest_df["consolidation_transcript_id"]
                                            == consol_transcript
                                        ]["consol_start"].values[0],
                                        int(row[2]),
                                        row[5],
                                        gene_id,
                                        consol_transcript,
                                    ]
                                ],
                                columns=consol_gene.columns,
                            ),
                        ],
                        ignore_index=True,
                    )
                elif num == len(bed_gene_data[transcript]):
                    consol_gene = pd.concat(
                        [
                            consol_gene,
                            pd.DataFrame(
                                [
                                    [
                                        row[0],
                                        int(row[1]) + 1,
                                        longest_df[
                                            longest_df["consolidation_transcript_id"]
                                            == consol_transcript
                                        ]["consol_end"].values[0],
                                        row[5],
                                        gene_id,
                                        consol_transcript,
                                    ]
                                ],
                                columns=consol_gene.columns,
                            ),
                        ],
                        ignore_index=True,
                    )
                else:
                    consol_gene = pd.concat(
                        [
                            consol_gene,
                            pd.DataFrame(
                                [
                                    [
                                        row[0],
                                        int(row[1]) + 1,
                                        int(row[2]),
                                        row[5],
                                        gene_id,
                                        consol_transcript,
                                    ]
                                ],
                                columns=consol_gene.columns,
                            ),
                        ],
                        ignore_index=True,
                    )
                    consol_gene = consol_gene.sort_values(
                        by=["transcript_id", "start"]
                    ).reset_index(drop=True)
    logger.debug(
        "Gene {} has {} transcript(s) before consolidation and {} transcript(s) after "
        "consolidation.",
        gene_id,
        len(key_gene),
        key_gene["consolidation_transcript_id"].nunique(),
    )
    return consol_gene, key_gene


def consolidate_transcripts(data, outdir, consol_prefix, consol_outfiles, genes, skip_interm):
    """
    If requested, consolidate junctions in input transcripts.
    """
    logger.info("Consolidation of transcript with identical junctions.")
    # Loop over genes
    consol_data = pd.DataFrame(columns=data.columns)
    if not skip_interm:
        consol_fhs = open_output_files(outdir, consol_outfiles)
        consol_fhs["key_fh"].write_text(",".join(consol_key_cols) + "\n")
    for gene in genes.groups:
        # Test for single transcript gene (WBGene00000003)
        # Test gene for multiple groups with consolidation (WBGene00001574)
        # Test monoexon gene (WBGene00000214)
        # Test monoexon with multiexon (WBGene00022161) -> added extra 2 monoexon transcripts
        # using gene_df =
        # pd.concat([gene_df,pd.DataFrame([["I",1779855,1781091,'-',"WBGene00022161","new_transcript_1"],["I",1781991,1782900,'-',"WBGene00022161","new_transcript_2"]],
        # columns=gene_df.columns)],ignore_index=True)
        pre_consol_jct = []
        gene_df = data[data["gene_id"] == gene]
        try:
            bed_gene_data = prep_bed_for_ea(gene_df)
        except ValueError as e:
            logger.error(e)
            continue
        transcripts = gene_df.groupby("transcript_id")
        # Loop over transcripts to get junctions
        for transcript in transcripts.groups:
            pre_consol_jct.extend(
                create_junction_catalog_str(gene, transcript, bed_gene_data[transcript])
            )
        pre_consol_jct_df = pd.DataFrame(pre_consol_jct, columns=jct_df_cols)
        # Consolidate 5'/3' variation
        consol_gene, key_gene = consolidate_junctions(
            bed_gene_data, pre_consol_jct_df, outdir, skip_interm, consol_prefix
        )
        consol_data = pd.concat([consol_data, consol_gene], ignore_index=True)
        if not skip_interm:
            write_output(key_gene, consol_fhs, "key_fh")
    # Output consolidated GTF
    if not skip_interm:
        write_gtf(consol_data, consol_fhs, "consol_gtf_fh")
    # Set data variable to new consolidated data
    data = consol_data.copy()
    del consol_data
    genes = data.groupby("gene_id")
    num_tx = data.groupby(["gene_id", "transcript_id"])
    logger.info(
        "After consolidation: {} genes and {} transcripts", len(genes), len(num_tx)
    )
    return data, genes

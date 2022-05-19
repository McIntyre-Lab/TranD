#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:33:04 2021

@author: adalena.nanni
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
from upsetplot import UpSet

#  Plot Functions


def plot_transcript_in_gene_split_pie(
    md_data, f1_odds, f2_odds, name1, name2, legendOut
):
    """
    Plot pie charts for the number of genes where the transcripts are the same,
            where one dataset has more transcripts that the onther, and where
            genes are exclusive (pie chart for all genes, single transcript in
            at least one dataset, multi-transcript in at least one dataset, and
            multi-transcript in both datasets
    """
    # Get counts
    cols = {"transcript_id": "num_transcript_in_gene_" + name1}
    geneF1only = (
        f1_odds.groupby("gene_id")["transcript_id"]
        .nunique()
        .reset_index()
        .rename(columns=cols)
    )

    geneF1only["num_transcript_in_gene_" + name2] = 0
    geneF1only["transcript_in_gene"] = name1 + "_only"
    cols = {"transcript_id": "num_transcript_in_gene_" + name2}
    geneF2only = (
        f2_odds.groupby("gene_id")["transcript_id"]
        .nunique()
        .reset_index()
        .rename(columns=cols)
    )
    geneF2only["num_transcript_in_gene_" + name1] = 0
    geneF2only["transcript_in_gene"] = name2 + "_only"

    # Make dictionary of names for plot and order
    geneCountNameDict = {
        "match": "Match",
        name1 + "_greater": name1 + " Greater",
        name2 + "_greater": name2 + " Greater",
        name1 + "_only": name1 + " Only",
        name2 + "_only": name2 + " Only",
    }

    geneCountOrderDict = {
        "match": 0,
        name1 + "_greater": 1,
        name2 + "_greater": 2,
        name1 + "_only": 3,
        name2 + "_only": 4,
    }

    # Plot pie chart
    # 1) all genes
    # 2) genes with 1 xcrpt for at least 1 dataset
    # 3) genes with >1 xcrpt for at least one dataset
    # 4) genes with > 1 xcrpt for both
    fig = plt.figure()
    for num in range(1, 5):
        ax = plt.subplot2grid((2, 2), (abs(num % 2 - 1), int(num / 3)), fig=fig)
        if num == 1:
            geneCount = pd.concat(
                [
                    md_data[["gene_id", "transcript_in_gene"]].drop_duplicates(),
                    geneF1only,
                    geneF2only,
                ],
                ignore_index=True,
                sort=False,
            )["transcript_in_gene"].value_counts(sort=False)
            title = "All Genes (n = {})".format(geneCount.sum())
            totalGenes = geneCount.sum()
            legendLabels = geneCount.reindex(list(geneCountOrderDict)).index.map(
                geneCountNameDict
            )
        elif num == 2:
            geneCount = pd.concat(
                [
                    md_data[
                        (md_data["num_transcript_in_gene_" + name1] == 1)
                        | (md_data["num_transcript_in_gene_" + name2] == 1)
                    ][["gene_id", "transcript_in_gene"]].drop_duplicates(),
                    geneF1only[geneF1only["num_transcript_in_gene_" + name1] == 1],
                    geneF2only[geneF2only["num_transcript_in_gene_" + name2] == 1],
                ],
                ignore_index=True,
                sort=False,
            )["transcript_in_gene"].value_counts(sort=False)
            title = "Genes with 1 Transcript\n(in at least one dataset, n = {})".format(
                geneCount.sum()
            )
        elif num == 3:
            geneCount = pd.concat(
                [
                    geneF1only[geneF1only["num_transcript_in_gene_" + name1] > 1],
                    geneF2only[geneF2only["num_transcript_in_gene_" + name2] > 1],
                ],
                ignore_index=True,
                sort=False,
            )["transcript_in_gene"].value_counts(sort=False)
            title = "Genes with >1 Transcript\n(exclusive to one dataset, n = {})".format(
                geneCount.sum()
            )
        else:
            geneCount = (
                md_data[
                    (md_data["num_transcript_in_gene_" + name1] > 1)
                    & (md_data["num_transcript_in_gene_" + name2] > 1)
                ][["gene_id", "transcript_in_gene"]]
                .drop_duplicates()["transcript_in_gene"]
                .value_counts(sort=False)
            )
            title = "Genes with >1 Transcript\n(in both datasets, n = {})".format(
                geneCount.sum()
            )
        geneCount.reindex(list(geneCountOrderDict)).plot(
            kind="pie",
            y="transcript_in_gene",
            labels=None,
            figsize=(12, 6),
            autopct=(lambda pct: get_pie_label(pct, geneCount)),
            colormap=ListedColormap(sns.color_palette("colorblind", 15).as_hex()),
            ax=ax,
        )
        ax.set_ylabel("")
        ax.set_title(title)
    plt.legend(bbox_to_anchor=(1.2, 1.2), loc="upper left", labels=legendLabels)
    plt.tight_layout()
    legendText = (
            "Percent of genes with equal number of transcripts between "
            "{} and {} (Match, blue), more transcripts in {} compared to {} "
            "({} Greater, dark orange), more transcripts in {} compared to {} "
            "({} Greater, gray), and genes exclusive to {} ({} Only, light "
            "orange) or {} ({} Only, pink). {} genes total.".format(
                name1,
                name2,
                name1,
                name2,
                name1,
                name2,
                name1,
                name2,
                name1,
                name1,
                name2,
                name2,
                int(totalGenes),
        )
    )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Shared and unique transcripts per gene \b0 \line "
            r"{}".format(
                legendText
            )
        )
        end_rtf(outFile)


def get_pie_label(pct, allvals):
    """
    Set pie label with percent and counts in ()
    """
    absolute = int(round(pct / 100.0 * np.sum(allvals)))
    return "{:.1f}%\n({:d})".format(pct, absolute)


def plot_gene_stack(md_data, name1, name2, legendOut, useProp=False):
    """
    Plot a stacked bar chart of the genes with reciprocal minimum pairs for all transcript pairs
    possible (Match: Reciprocal Pairs or Greater: Reciprocal Pairs), only some of the transcript
    pairs (Match: Partial Reciprocal Pairs or Greater: Partial Reciprocal Pairs), or none of the
    transcript pars (Match: No Reciprocal Pairs or Greater: No Reciprocal Pairs) using gene counts
    or proportions of Match or Greater genes
    """
    genePairDF = (
        md_data.groupby(["transcript_in_gene", "recip_min_pair_in_gene"])["gene_id"]
        .nunique()
        .reset_index()
        .pivot_table(
            index=["transcript_in_gene"],
            columns="recip_min_pair_in_gene",
            values="gene_id",
        )
    )
    for col in ["reciprocal_pairs", "partial_reciprocal_pairs", "no_reciprocal_pairs"]:
        if col not in genePairDF.columns:
            genePairDF[col] = 0
    genePairDF = genePairDF.fillna(0)
    genePairDF = genePairDF.rename(
        columns={
            "reciprocal_pairs": "Reciprocal Pairs",
            "partial_reciprocal_pairs": "Partial Reciprocal Pairs",
            "no_reciprocal_pairs": "No Reciprocal Pairs",
        },
        index={
            name1 + "_greater": name1 + " Greater",
            name2 + "_greater": name2 + " Greater",
            "match": "Match",
        },
    )
    genePairPallete = ["#FAC748", "#1D2F6F", "#8390FA"]
    if useProp:
        legendText = (
            "The stacked bar chart is scaled to represent proportions of "
            "genes within each category. Columns are: genes with equal numbers "
            "of transcripts in {} and {} (Match; n = {}), genes with {} "
            "containing more transcripts than {} ({} Greater; n = {}), and "
            "genes with {} containing more transcripts than {} ({} Greater, "
            "n = {}). Genes with Reciprocal Pairs for all transcripts (Match) "
            "or as a complete subset are colored blue. Genes with at least one "
            "but not all transcript pairs are reciprocal minimum matches "
            "(Partial Reciprocal Pairs) are colored yellow. Genes with no "
            "pairs that are reciprocal minimum matches (No Reciprocal Pairs) "
            "are colored purple.".format(
                md_data["gene_id"].nunique(),
                name1,
                name2,
                genePairDF[genePairDF.index == "Match"].sum().sum(),
                name1,
                name2,
                name1,
                genePairDF[genePairDF.index == name1 + " Greater"].sum().sum(),
                name2,
                name1,
                name2,
                genePairDF[genePairDF.index == name2 + " Greater"].sum().sum(),
            )
        )
        genePairDF = genePairDF.div(genePairDF.sum(axis=1), axis=0)
        title = "Proportion of Genes"
    else:
        legendText = (
            "The stacked bar chart is the count of the number of genes in "
            "each category (N = {}). Columns are: genes with equal numbers of "
            "transcripts in {} and {} (Match; n = {}), genes with {} "
            "containing more transcripts than {} ({} Greater; n = {}), and "
            "genes with {} containing more transcripts than {} ({} Greater, "
            "n = {}). Genes with Reciprocal Pairs for all transcripts (Match) "
            "or as a complete subset are colored blue. Genes with at least one "
            "but not all transcript pairs are reciprocal minimum matches "
            "(Partial Reciprocal Pairs) are colored yellow. Genew with no "
            "pairs that are reciprocal minimum matches (No Reciprocal Pairs) "
            "are colored purple.".format(
                md_data["gene_id"].nunique(),
                name1,
                name2,
                genePairDF[genePairDF.index == "Match"].sum().sum(),
                name1,
                name2,
                name1,
                genePairDF[genePairDF.index == name1 + " Greater"].sum().sum(),
                name2,
                name1,
                name2,
                genePairDF[genePairDF.index == name2 + " Greater"].sum().sum(),
            )
        )
        title = "Number of Genes"
    genePairDF.plot(
        kind="bar",
        figsize=(8, 6),
        stacked=True,
        rot=45,
        colormap=ListedColormap(genePairPallete),
    )
    plt.ylabel(title)
    plt.xlabel("")
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.tight_layout()
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Classification of transcript pairs \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_transcript_in_gene_upset(md_data, f1_odds, f2_odds, name1, name2, legendOut):
    """
    UpSet plot for number of transcripts in genes
    """
    geneF1only = (
        f1_odds.groupby("gene_id")["transcript_id"]
        .nunique()
        .reset_index()
        .rename(columns={"transcript_id": "num_transcript_in_gene_" + name1})
    )
    geneF1only["num_transcript_in_gene_" + name2] = 0
    geneF1only["transcript_in_gene"] = name1 + "_only"
    geneF2only = (
        f2_odds.groupby("gene_id")["transcript_id"]
        .nunique()
        .reset_index()
        .rename(columns={"transcript_id": "num_transcript_in_gene_" + name2})
    )
    geneF2only["num_transcript_in_gene_" + name1] = 0
    geneF2only["transcript_in_gene"] = name2 + "_only"
    geneAll = pd.concat(
        [
            md_data[
                [
                    "gene_id",
                    "num_transcript_in_gene_" + name1,
                    "num_transcript_in_gene_" + name2,
                    "transcript_in_gene",
                ]
            ].drop_duplicates(),
            geneF1only,
            geneF2only,
        ],
        ignore_index=True,
        sort=False,
    )
    colList = []
    # Flag when genes have 1-4 transcripts
    for num in range(1, 5):
        geneAll[str(num) + " Transcript(s) " + name1] = np.where(
            geneAll["num_transcript_in_gene_" + name1] == num, True, False
        )
        colList.append(str(num) + " Transcript(s) " + name1)
        geneAll[str(num) + " Transcript(s) " + name2] = np.where(
            geneAll["num_transcript_in_gene_" + name2] == num, True, False
        )
        colList.append(str(num) + " Transcript(s) " + name2)
    # Flag when genes have 5+ transcripts
    geneAll["5+ Transcript(s) " + name1] = np.where(
        geneAll["num_transcript_in_gene_" + name1] >= 5, True, False
    )
    colList.append("5+ Transcript(s) " + name1)
    geneAll["5+ Transcript(s) " + name2] = np.where(
        geneAll["num_transcript_in_gene_" + name2] >= 5, True, False
    )
    colList.append("5+ Transcript(s) " + name2)
    legendText = (
        "Rows indicate the number of transcripts in {} and {}. Columns are "
        "separated by unique combinations of the number of transcripts in "
        "{} and {}. When a gene is not present in one of the GTF files, only "
        "a single dot will be indicated. The number of genes (N = {}) for "
        "each nonoverlapping category are displayed in a histogram.".format(
            name1,
            name2,
            name1,
            name2,
            sum(
                [
                    get_value_count(geneAll, "transcript_in_gene", name1 + "_only"),
                    get_value_count(geneAll, "transcript_in_gene", name2 + "_only"),
                    get_value_count(geneAll, "transcript_in_gene", "match"),
                    get_value_count(geneAll, "transcript_in_gene", name1 + "_greater"),
                    get_value_count(geneAll, "transcript_in_gene", name2 + "_greater"),                    
                ]
            ),
        )
    )
    plot_upset(
        geneAll.set_index(colList), "Number of Genes with Each Number of Transcripts"
    )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Number of transcripts per gene comparison UpSet plot"
            r"\b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_min_pair_AS_upset_nt_box(md_data, name1, name2, legendOut, reciprocal=True, pairs=None):
    """
    For all reciprocal minimum pair transcripts or minimum pair of extra transcripts:
            Upset plot for the number of transcript pairs with each kind of
            alternative splicing (AS) and the distribution of nt differences
            between the pairs plotted
    """
    if reciprocal:
        # Plot only reciprocal minimum pairs
        minPairAS = md_data[md_data["flag_recip_min_match"] == 1].copy()
    else:
        if pairs == None:
            # Plot only minimum pairs of extras, or those without a recip min pair
            minPairAS = md_data[
                (md_data["flag_recip_min_match"] != 1)
                & (
                    (md_data["flag_min_match_" + name1] == 1)
                    | (md_data["flag_min_match_" + name2] == 1)
                )
            ].copy()
        elif pairs == "all":
            minPairAS = md_data[
                    (md_data["flag_min_match_" + name1] == 1)
                    | (md_data["flag_min_match_" + name2] == 1)
            ].copy()
        elif pairs == "first":
            minPairAS = md_data[
                    (md_data["flag_min_match_" + name1] == 1)
            ].copy()
        elif pairs == "second":
            minPairAS = md_data[
                    (md_data["flag_min_match_" + name2] == 1)
            ].copy()
    # Check if there are any pairs to plot
    if len(minPairAS) == 0:
        return
    minPairAS = minPairAS[
        [
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
    ]
    # Set No Shared NT genes to have np.nan nucleotide difference
    minPairAS["num_nt_diff"] = np.where(
        minPairAS["flag_no_shared_nt"] == 1,
        np.nan,
        minPairAS["num_nt_diff"]
    )
    minPairAS = minPairAS.rename(
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
    AScols = [
        "5' Variation",
        "3' Variation",
        "Alt. Donor/Acceptor",
        "Alt. Exon",
        "Intron Retention",
        "No Shared NT",
    ]
    minPairAS[AScols] = minPairAS[AScols].astype(int)
    if reciprocal:
        plot_upset(
            minPairAS.set_index(AScols),
            "Number of Reciprocal Minimum Pairs with AS Categories",
            ["# NT Different", "Proportion\nNT Different"],
        )
        legendText = (
            "Number of reciprocal minimum pairs with the specified types of alternative "
            "splicing between {} and {} indicated by the black dots below the histogram "
            "of pair counts (n = {} total pairs). Pairs with no shared nucleotides are "
            "in the same genes but with nonoverlapping coordinates. Box plots of the "
            "number (blue) and proportion (orange) of nucleotide (NT) differences between "
            "the pairs represented in the histogram.".format(
                name1,
                name2,
                len(minPairAS),
            )
        )
    else:
        if pairs == None:
            plot_upset(
                minPairAS.set_index(AScols),
                "Number of Extra Minimum Pairs with AS Categories",
                ["# NT Different", "Proportion\nNT Different"],
            )
            legendText = (
                "Number of extra minimum pairs, or pairs without reciprocal minimums, "
                "with the specified types of alternative splicing between {} and {} "
                "indicated by the black dots below the histogram of pair counts "
                "(n = {} total pairs). Pairs with no shared nucleotides are "
                "in the same genes but with nonoverlapping coordinates. Box plots of the "
                "number (blue) and proportion (orange) of nucleotide (NT) differences between "
                "the pairs represented in the histogram.".format(
                    name1,
                    name2,
                    len(minPairAS),
                )
            )
        elif pairs == "all":
            plot_upset(
                minPairAS.set_index(AScols),
                "Number of Minimum Pairs in either {} or {} with AS Categories".format(name1, name2),
                ["# NT Different", "Proportion\nNT Different"],
            )
            legendText = (
                "Number of minimum pairs in either {} or {}  "
                "with the specified types of alternative splicing between {} and {} "
                "indicated by the black dots below the histogram of pair counts "
                "(n = {} total pairs). Pairs with no shared nucleotides are "
                "in the same genes but with nonoverlapping coordinates. Box plots of the "
                "number (blue) and proportion (orange) of nucleotide (NT) differences between "
                "the pairs represented in the histogram.".format(
                    name1,
                    name2,
                    name1,
                    name2,
                    len(minPairAS),
                )
            )
        elif pairs == "first":
            plot_upset(
                minPairAS.set_index(AScols),
                "Number of Minimum Pairs in {} with AS Categories".format(name1),
                ["# NT Different", "Proportion\nNT Different"],
            )
            legendText = (
                "Number of minimum pairs in {}  "
                "with the specified types of alternative splicing between {} and {} "
                "indicated by the black dots below the histogram of pair counts "
                "(n = {} total pairs). Pairs with no shared nucleotides are "
                "in the same genes but with nonoverlapping coordinates. Box plots of the "
                "number (blue) and proportion (orange) of nucleotide (NT) differences between "
                "the pairs represented in the histogram.".format(
                    name1,
                    name1,
                    name2,
                    len(minPairAS),
                )
            )
        elif pairs == "second":
            plot_upset(
                minPairAS.set_index(AScols),
                "Number of Minimum Pairs in {} with AS Categories".format(name2),
                ["# NT Different", "Proportion\nNT Different"],
            )
            legendText = (
                "Number of minimum pairs in {}  "
                "with the specified types of alternative splicing between {} and {} "
                "indicated by the black dots below the histogram of pair counts "
                "(n = {} total pairs). Pairs with no shared nucleotides are "
                "in the same genes but with nonoverlapping coordinates. Box plots of the "
                "number (blue) and proportion (orange) of nucleotide (NT) differences between "
                "the pairs represented in the histogram.".format(
                    name2,
                    name1,
                    name2,
                    len(minPairAS),
                )
            )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        if reciprocal:
            outFile.write(
                r"\b Figure. Reciprocal minimum pair alternative splicing \b0"
                r" \line {}".format(
                    legendText
                )
            )
        else:
            outFile.write(
                r"\b Figure. Extra minimum pair alternative splicing \b0 \line "
                r"{}".format(
                    legendText
                )
            )
        end_rtf(outFile)


def plot_gene_recip_min_AS_upset(md_data, name1, name2, legendOut):
    """
    Upset plot for number of genes with each kind of AS when comparing the
            reciprocally min match pairs of the two datasets with added box plots
            of the number of min match pairs present in the gene and the number
            of transcripts in each dataset
    """
    recipMinPairAS = md_data[md_data["flag_recip_min_match"] == 1][
        [
            "gene_id",
            "flag_alt_exon_recip_min_match",
            "flag_alt_donor_acceptor_recip_min_match",
            "flag_IR_recip_min_match",
            "flag_5_variation_recip_min_match",
            "flag_3_variation_recip_min_match",
            "flag_no_shared_nt_recip_min_match",
            "num_recip_min_match_in_gene",
            "num_transcript_in_gene_" + name1,
            "num_transcript_in_gene_" + name2,
        ]
    ].copy()
    geneRecipMatchAS = (
        recipMinPairAS.groupby("gene_id")
        .agg(
            {
                "flag_alt_exon_recip_min_match": "max",
                "flag_alt_donor_acceptor_recip_min_match": "max",
                "flag_IR_recip_min_match": "max",
                "flag_5_variation_recip_min_match": "max",
                "flag_3_variation_recip_min_match": "max",
                "flag_no_shared_nt_recip_min_match": "max",
                "num_recip_min_match_in_gene": "max",
                "num_transcript_in_gene_" + name1: "max",
                "num_transcript_in_gene_" + name2: "max",
            }
        )
        .reset_index()
    )
    geneFlagCols = [
        "flag_alt_exon_recip_min_match",
        "flag_alt_donor_acceptor_recip_min_match",
        "flag_IR_recip_min_match",
        "flag_5_variation_recip_min_match",
        "flag_3_variation_recip_min_match",
        "flag_no_shared_nt_recip_min_match",
    ]
    geneRecipMatchAS[geneFlagCols] = geneRecipMatchAS[geneFlagCols].astype(bool)
    geneRecipMatchAS = geneRecipMatchAS.rename(
        columns={
            "flag_alt_exon_recip_min_match": "Alt. Exon",
            "flag_alt_donor_acceptor_recip_min_match": "Alt. Donor/Acceptor",
            "flag_IR_recip_min_match": "Intron Retention",
            "flag_5_variation_recip_min_match": "5' Variation",
            "flag_3_variation_recip_min_match": "3' Variation",
            "flag_no_shared_nt_recip_min_match": "No Shared NT",
            "num_recip_min_match_in_gene": "# Recip. Min.\nTranscript Pairs",
            "num_transcript_in_gene_" + name1: "# Transcripts\nin " + name1,
            "num_transcript_in_gene_" + name2: "# Transcripts\nin " + name2,
        }
    )
    AScols = [
        "5' Variation",
        "3' Variation",
        "Alt. Donor/Acceptor",
        "Alt. Exon",
        "Intron Retention",
        "No Shared NT",
    ]
    plot_upset(
        geneRecipMatchAS.set_index(AScols),
        "",
        [
            "# Recip. Min.\nTranscript Pairs",
            "# Transcripts\nin " + name1,
            "# Transcripts\nin " + name2,
        ],
    )
    legendText = (
        "Number of genes with the specified types of alternative splicing in only "
        "reciprocal minimum pairs between {} and {} indicated by the black dots "
        "below the histogram of gene counts (n = {} genes with {} reciprocal "
        "minimum pairs). Box plots represent the number of reciprocal minimum "
        "pairs (blue) and the number of transcripts in {} (orange) and {} "
        "(green). Genes with \"No Shared NT\" have a pair of "
        "transcripts with nonoverlapping coordinates.".format(
            name1,
            name2,
            len(geneRecipMatchAS),
            len(md_data[md_data["flag_recip_min_match"] == 1]),
            name1,
            name2,
        )
    )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Observed alternative splicing in genes among "
            r"reciprocal minimum pairs of transcripts \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_gene_AS_upset(td_data, legendOut):
    """
    Upset plot for number of genes with each kind of AS when comparing all
            transcript pairs within each gene of a single transcriptome with added box plot
            of the number transcripts per gene
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
        ]
    ].copy()
    # Get AS flags at the gene level
    geneAS = (
        pairAS.groupby("gene_id")
        .agg(
            {
                "flag_alt_exon": "max",
                "flag_alt_donor_acceptor": "max",
                "flag_IR": "max",
                "flag_5_variation": "max",
                "flag_3_variation": "max",
                "flag_no_shared_nt": "max",
            }
        )
        .reset_index()
    )
    # Get the number of transcripts per gene by getting the unique list
    #   of transcript_1 and transcript_2
    transcriptPerGene = (
        pd.concat(
            [
                pairAS[["gene_id", "transcript_1"]].rename(
                    columns={"transcript_1": "transcript_id"}
                ),
                pairAS[["gene_id", "transcript_2"]].rename(
                    columns={"transcript_2": "transcript_id"}
                ),
            ]
        )
        .groupby("gene_id")["transcript_id"]
        .nunique()
        .reset_index()
        .rename(columns={"transcript_id": "num_transcript_per_gene"})
    )
    # Merge transcripts per gene with AS flags
    mergeASxcrptPerGene = pd.merge(
        geneAS, transcriptPerGene, how="outer", on="gene_id", indicator="merge_check"
    )
    # Set AS flags to boolean values and rename
    geneFlagCols = [
        "gene_id",
        "flag_alt_exon",
        "flag_alt_donor_acceptor",
        "flag_IR",
        "flag_5_variation",
        "flag_3_variation",
        "flag_no_shared_nt",
    ]
    mergeASxcrptPerGene[geneFlagCols] = mergeASxcrptPerGene[geneFlagCols].astype(bool)
    mergeASxcrptPerGene = mergeASxcrptPerGene.rename(
        columns={
            "flag_alt_exon": "Alt. Exon",
            "flag_alt_donor_acceptor": "Alt. Donor/Acceptor",
            "flag_IR": "Intron Retention",
            "flag_5_variation": "5' Variation",
            "flag_3_variation": "3' Variation",
            "flag_no_shared_nt": "No Shared NT",
            "num_transcript_per_gene": "# Transcripts\nPer Gene",
            "num_nt_diff": "Avg #\nNT Different",
            "prop_nt_diff": "Avg Proportion\nNT Different",
        }
    )
    AScols = [
        "5' Variation",
        "3' Variation",
        "Alt. Donor/Acceptor",
        "Alt. Exon",
        "Intron Retention",
        "No Shared NT",
    ]
    plot_upset(
        mergeASxcrptPerGene.set_index(AScols),
        "",
        [
            "# Transcripts\nPer Gene",
        ],
    )
    legendText = (
        "Number of genes with the specified types of alternative splicing indicated "
        "by the black dots below the histogram of gene counts (n = {} multi-transcript genes). "
        "The box plot represent the number of transcripts per gene (blue). "
        "Genes with \"No Shared NT\" have a pair of transcripts with "
        "nonoverlapping coordinates.".format(
            len(mergeASxcrptPerGene)
        )
    )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Alternative splicing in genes \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_pair_AS_upset_nt_box(td_data, legendOut):
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
    AScols = [
        "5' Variation",
        "3' Variation",
        "Alt. Donor/Acceptor",
        "Alt. Exon",
        "Intron Retention",
        "No Shared NT",
    ]
    plot_upset(
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
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Alternative splicing between pairs of transcripts \b0"
            r" \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_upset(df, title, boxCols=None):
    """
    Plot an UpSet plot given a properly formatted dataframe (multi-indexed with boolean values)
    """
    upset = UpSet(
        df,
        subset_size="count",
        show_counts=True,
        sort_by="degree",
        sort_categories_by=None,
    )
    if boxCols is not None:
        if type(boxCols) != list:
            boxCols = list(boxCols)
        for col in boxCols:
            upset.add_catplot(
                value=col,
                kind="box",
                elements=4,
                showfliers=False,
                color=sns.color_palette("colorblind", 15).as_hex()[boxCols.index(col)],
            )
    upset.plot()
    plt.subplots_adjust(right=1.00001)
    plt.suptitle(title)


def plot_gene_avg_nt_diff_pairs(md_data, name1, name2, legendOut, zoomMean=False):
    """
    Plot average number of nt different in recip. min match pairs against
            avg number of nt different in min pairs that are extra
            (either in the greater sets of in set that do not match) with the
            option to zoom in to the  mean value of non-zero nucleotide differences
            (max of the two means from recip. min pairs and extras)
    """
    minPairNT = md_data[
        [
            "gene_id",
            "transcript_in_gene",
            "recip_min_pair_in_gene",
            "flag_recip_min_match",
            "flag_min_match_" + name1,
            "flag_min_match_" + name2,
            "num_nt_diff",
            "prop_nt_diff",
        ]
    ].copy()
    # Ensure number of nt different are int and float values
    minPairNT["num_nt_diff"] = minPairNT["num_nt_diff"].astype(int)
    minPairNT["prop_nt_diff"] = minPairNT["prop_nt_diff"].astype(float)
    minPairNT["num_nt_diff_recip_min"] = np.where(
        minPairNT["flag_recip_min_match"] == 1, minPairNT["num_nt_diff"], 0
    )
    minPairNT["num_nt_diff_extra"] = np.where(
        (minPairNT["flag_recip_min_match"] == 0)
        & (
            minPairNT["flag_min_match_" + name1] + minPairNT["flag_min_match_" + name2]
            == 1
        ),
        minPairNT["num_nt_diff"],
        0,
    )
    minPairNTGene = (
        minPairNT.groupby("gene_id")
        .agg(
            {
                "transcript_in_gene": "first",
                "recip_min_pair_in_gene": "first",
                "num_nt_diff_recip_min": "mean",
                "num_nt_diff_extra": "mean",
            }
        )
        .reset_index()
    )
    recipConditions = [
        (minPairNTGene["transcript_in_gene"] == "match")
        & (minPairNTGene["recip_min_pair_in_gene"] == "reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == "match")
        & (minPairNTGene["recip_min_pair_in_gene"] == "partial_reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == "match")
        & (minPairNTGene["recip_min_pair_in_gene"] == "no_reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == name1 + "_greater")
        & (minPairNTGene["recip_min_pair_in_gene"] == "reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == name1 + "_greater")
        & (minPairNTGene["recip_min_pair_in_gene"] == "partial_reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == name1 + "_greater")
        & (minPairNTGene["recip_min_pair_in_gene"] == "no_reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == name2 + "_greater")
        & (minPairNTGene["recip_min_pair_in_gene"] == "reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == name2 + "_greater")
        & (minPairNTGene["recip_min_pair_in_gene"] == "partial_reciprocal_pairs"),
        (minPairNTGene["transcript_in_gene"] == name2 + "_greater")
        & (minPairNTGene["recip_min_pair_in_gene"] == "no_reciprocal_pairs"),
    ]
    recipChoices = [
        "Match:Reciprocal Pairs",
        "Match:Partial Reciprocal Pairs",
        "Match:No Reciprocal Pairs",
        name1 + " Greater:Reciprocal Pairs",
        name1 + " Greater:Partial Reciprocal Pairs",
        name1 + " Greater:No Reciprocal Pairs",
        name2 + " Greater:Reciprocal Pairs",
        name2 + " Greater:Partial Reciprocal Pairs",
        name2 + " Greater:No Reciprocal Pairs",
    ]
    minPairNTGene["Category"] = np.select(recipConditions, recipChoices, "NC")
    minPairGeneOrder = {
        "Match:Reciprocal Pairs": 0,
        "Match:Partial Reciprocal Pairs": 1,
        "Match:No Reciprocal Pairs": 2,
        name1 + " Greater:Reciprocal Pairs": 3,
        name1 + " Greater:Partial Reciprocal Pairs": 4,
        name1 + " Greater:No Reciprocal Pairs": 5,
        name2 + " Greater:Reciprocal Pairs": 6,
        name2 + " Greater:Partial Reciprocal Pairs": 7,
        name2 + " Greater:No Reciprocal Pairs": 8,
    }
    minPairNTGene["Category"] = pd.Categorical(
        minPairNTGene["Category"],
        categories=sorted(minPairGeneOrder, key=minPairGeneOrder.get),
        ordered=True,
    )
    minPairNTGene = minPairNTGene.sort_values("Category")
    colorPalleteRecip = {
        "Match:Reciprocal Pairs": "#0173b2",
        "Match:Partial Reciprocal Pairs": "#56b4e9",
        "Match:No Reciprocal Pairs": "#029e73",
        name1 + " Greater:Reciprocal Pairs": "#d55e00",
        name1 + " Greater:Partial Reciprocal Pairs": "#de8f05",
        name1 + " Greater:No Reciprocal Pairs": "#ece133",
        name2 + " Greater:Reciprocal Pairs": "#949494",
        name2 + " Greater:Partial Reciprocal Pairs": "#fbafe4",
        name2 + " Greater:No Reciprocal Pairs": "#ca9161",
    }
    if not zoomMean:
        limit = (
            max(
                minPairNTGene[
                        ["num_nt_diff_recip_min", "num_nt_diff_extra"]
                        ].max().max(),
                1
            )
        )
    else:
        try:
            limit = round(
                np.nanmax(
                    [
                        minPairNTGene[minPairNTGene["num_nt_diff_recip_min"] > 0][
                            "num_nt_diff_recip_min"
                        ].mean(),
                        minPairNTGene[minPairNTGene["num_nt_diff_extra"] > 0][
                            "num_nt_diff_extra"
                        ].mean(),
                    ]
                )
            )
        except ValueError:
            limit = 1
    plt.figure(figsize=(8.45, 5))
    sns.scatterplot(
        data=minPairNTGene,
        x="num_nt_diff_recip_min",
        y="num_nt_diff_extra",
        hue="Category",
        palette=colorPalleteRecip,
    )
    plt.xlim(0, limit)
    plt.ylim(0, limit)
    plt.xlabel("Avg. # NT Different in Recip. Min. Pairs")
    plt.ylabel("Avg. # NT Different in Min. Pairs of Extras")
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.plot((0, limit), (0, limit), color="r")
    legendText = (
        "Each point represents a gene plotted by the average number of nucleotides "
        "different between reciprocal minimum pairs (X-axis) and minimum pairs "
        "of extras without a reciprocal minimum pair (Y-axis). Genes are colored "
        "by categorizations of the comparison between {} and {}.".format(
            name1,
            name2,
        )
    )
    plt.tight_layout()
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Number of nucleotides different between reciprocal "
            r"and nonreciprocal minimum pairs \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_complexity_box(geneDF, transcriptDF, outdir, legendOut):
    """
    Plot box plots of complexity measures:
        1) Transcripts per gene
        2) Unique exons per gene
        3) Exons per transcript
    """
    #    fig = plt.figure(figsize=(6,12))
    fig = plt.figure(figsize=(12, 5))
    #    axTop = plt.subplot2grid((3,1),(0,0),fig=fig)
    axTop = plt.subplot2grid((1, 3), (0, 0), fig=fig)
    axTop.set_xlabel("")
    #    axTop.set_ylabel("Number of Transcripts Per Gene")
    axTop.set_title("Number of Transcripts Per Gene")
    axTop.set_ylabel("")
    #    axMid = plt.subplot2grid((3,1),(1,0),fig=fig)
    axMid = plt.subplot2grid((1, 3), (0, 1), fig=fig)
    axMid.set_xlabel("")
    #    axMid.set_ylabel("Number of Unique Exons Per Gene")
    axMid.set_title("Number of Unique Exons Per Gene")
    axMid.set_ylabel("")
    #    axBottom = plt.subplot2grid((3,1),(2,0),fig=fig)
    axBottom = plt.subplot2grid((1, 3), (0, 2), fig=fig)
    axBottom.set_xlabel("")
    #    axBottom.set_ylabel("Number of Exons Per Transcript")
    axBottom.set_title("Number of Exons Per Transcript")
    axBottom.set_ylabel("")
    bplotTop = geneDF.set_index("gene_id")[["num_transcript"]].boxplot(
        ax=axTop, notch=True, patch_artist=True, return_type="dict"
    )
    axTop.set_xticks([])
    bplotMid = geneDF.set_index("gene_id")[["num_uniq_exon"]].boxplot(
        ax=axMid, notch=True, patch_artist=True, return_type="dict"
    )
    axMid.set_xticks([])
    bplotBottom = transcriptDF.set_index("transcript_id")[["num_exon"]].boxplot(
        ax=axBottom, notch=True, patch_artist=True, return_type="dict"
    )
    axBottom.set_xticks([])
    for bplot in (bplotTop, bplotMid, bplotBottom):
        for patch in bplot["boxes"]:
            patch.set_facecolor("turquoise")
            patch.set_edgecolor("black")
        for median in bplot["medians"]:
            median.set(color="red")
        for whisker in bplot["whiskers"]:
            whisker.set(color="black")
    legendText = (
        "Distributions of transcriptome complexity measures including transcripts "
        "per gene (median = {}, mean = {}), unique exons per gene "
        "(median = {}, mean = {}), and exons per transcript "
        "(median = {}, mean = {}). Total number of genes = {}. "
        "Total number of transcripts = {}.".format(
            geneDF["num_transcript"].median(),
            geneDF["num_transcript"].mean(),
            geneDF["num_uniq_exon"].median(),
            geneDF["num_uniq_exon"].mean(),
            transcriptDF["num_exon"].median(),
            transcriptDF["num_exon"].mean(),
            len(geneDF),
            int(geneDF["num_transcript"].sum()),
        )
    )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Transcriptome complexity \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_complexity_violin(geneDF, transcriptDF, outdir, legendOut):
    """
    Plot violin plots of complexity measures:
        1) Transcripts per gene
        2) Unique exons per gene
        3) Exons per transcript
    """
    sns.set(style="whitegrid")
    fig = plt.figure(figsize=(12, 5))
    axLeft = plt.subplot2grid((1, 3), (0, 0), fig=fig)
    axLeft.set_title("Number of Transcripts Per Gene")
    axMid = plt.subplot2grid((1, 3), (0, 1), fig=fig)
    axMid.set_title("Number of Unique Exons Per Gene")
    axRight = plt.subplot2grid((1, 3), (0, 2), fig=fig)
    axRight.set_title("Number of Exons Per Transcript")
    sns.violinplot(y="num_transcript", data=geneDF.set_index("gene_id"), ax=axLeft)
    sns.violinplot(y="num_uniq_exon", data=geneDF.set_index("gene_id"), ax=axMid)
    sns.violinplot(
        y="num_exon", data=transcriptDF.set_index("transcript_id"), ax=axRight
    )
    axLeft.set_xlabel("")
    axLeft.set_ylabel("")
    axMid.set_xlabel("")
    axMid.set_ylabel("")
    axRight.set_xlabel("")
    axRight.set_ylabel("")
    plt.tight_layout()
    legendText = (
        "Distributions of transcriptome complexity measures including transcripts "
        "per gene (median = {}, mean = {}), unique exons per gene (median = {}, mean = {}), "
        "and exons per transcript (median = {}, mean = {}). "
        "Total genes described = {}.".format(
            geneDF["num_transcript"].median(),
            geneDF["num_transcript"].mean(),
            geneDF["num_uniq_exon"].median(),
            geneDF["num_uniq_exon"].mean(),
            transcriptDF["num_exon"].median(),
            transcriptDF["num_exon"].mean(),
            len(geneDF),
        )
    )
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Transcriptome complexity \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def plot_gene_prop_nt_variablility(ef_data, legendOut, multitranscript=False):
    """
    Plot distribution of proportion of nucleotide variability across all genes or
    only multitranscript genes
    """
    # Get lengths of exon fragments
    ef_data["num_nt_varying"] = np.where(
        ef_data["ea_annotation_frequency"] != "constitutive",
        ef_data["ef_end"] - ef_data["ef_start"],
        0,
    )
    ef_data["num_nt_total"] = ef_data["ef_end"] - ef_data["ef_start"]
    # For each gene, get the number of varying nucleotides (not constitutive)
    ef_gene_data = (
        ef_data.groupby("gene_id")[["num_nt_varying", "num_nt_total"]]
        .sum()
        .reset_index()
    )
    # Get propotion of varying nucleotides
    ef_gene_data["prop_varying_nt"] = (
        ef_gene_data["num_nt_varying"] / ef_gene_data["num_nt_total"]
    )

    # Check if outputting all genes or only multitranscript genes
    if multitranscript:
        ef_data_split = split_column_by_sep(
            ef_data,
            col_name="ef_transcript_ids",
            sep="|"
        )
        xcrpt_per_gene = (
            ef_data_split.groupby("gene_id")["ef_transcript_ids"]
            .nunique()
            .reset_index()
        )
        multi_xcrpt_gene = xcrpt_per_gene[xcrpt_per_gene["ef_transcript_ids"] > 1]
        ef_gene_data = ef_gene_data[
            ef_gene_data["gene_id"].isin(multi_xcrpt_gene["gene_id"])
        ]
        legendText = (
            "Distribution of the proportion of variable "
            "nucleotides across multi-transcript genes (n = {}).".format(
                len(ef_gene_data)
            )
        )
    else:
        legendText = (
            "Distribution of the proportion of variable "
            "nucleotides across all genes (n = {}). "
            "Single transcript genes are assigned a proportion of 1.".format(
                len(ef_gene_data)
            )
        )
    # Plot density
    sns.distplot(
            ef_gene_data["prop_varying_nt"],
            hist=False
    )
    plt.xlim(0, 1)
    plt.ylabel("Density")
    plt.xlabel("Proportion of Variable Nucleotides")
    plt.tight_layout()
    with open(legendOut, "w") as outFile:
        start_rtf(outFile)
        outFile.write(
            r"\b Figure. Transcriptome complexity \b0 \line {}".format(
                legendText
            )
        )
        end_rtf(outFile)


def split_column_by_sep(df, col_name=None, sep=None, sort_list=None):
    """
    Split variable by some character like '|' or ',' and
    keep all other values the same
    """
    if col_name is None:
        col_name = "transcript_id"
    if sep is None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del tempDF[col_name]
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list is not None:
        splitDF = splitDF.sort_values(by=sort_list)
    del (tempDF, splitList)
    return splitDF


def get_value_count(df, col, value):
    """
    Do value_count function of column col on dataframe df
    Get count of value (an element in col), if value is not present then return 0
    """
    df2 = df[col].value_counts()
    if value in df2.index:
        value = df2[value]
    else:
        value = 0
    return value


def start_rtf(outFile):
    """
    Open RTF file
    """
    outFile.write(
        r"{\rtf1\ansi\ansicpg1252\deff0\deflang1033{\fonttbl{\f0\fswiss\fcharset0 Arial;}}"
    )

def get_citation():
    return (
        r" \line \line 1. Nanni, A., Titus-McQuillan, J., Moskalenko, O., "
        r"Pardo-Palacios, F., Liu, Z., Conesa, A., Rogers, R. L., & McIntyre, "
        "L. M. (2021). The evolution of splicing: transcriptome complexity "
        "and transcript distances implemented in TranD. \i bioRxiv\i0, "
        "2021.2009.2028.462251. https://doi.org/10.1101/2021.09.28.462251. "
        r"https://github.com/McIntyre-Lab/TranD."
    )

def end_rtf(outFile):
    """
    Close RTF file
    """
    outFile.write(
        r" Transcriptome analyses performed by TranD [1]. {}"
        r" \line \line \line \i Disclaimer:  While automated captions of "
        r"TranD have been carefully constructed, users are advised to verify "
        r"caption contents before use and report any errors to the TranD "
        r"github.\i0 ".format(
            get_citation()
        )
    )
    outFile.write(r"}\n\x00")

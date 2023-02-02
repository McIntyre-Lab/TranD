#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: adalena
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Functions for TranD plots
import seaborn as sns
from upsetplot import UpSet


testdf = pd.read_csv(
    "C:/Users/knife/OneDrive - University of Florida/McIntyre Lab/MyTranD/mysrc/my_test_csv.csv")

plot_pair_AS_upset_nt_box(testdf, "C:/Users/knife/Desktop", True)


def plot_pair_AS_upset_nt_box(td_data, legendOut, drop_5_3_no_shared_nt=False):
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
    if drop_5_3_no_shared_nt:
        AScolsSubset = [
            "Alt. Donor/Acceptor",
            "Alt. Exon",
            "Intron Retention",
        ]
        plot_upset(
            pairAS.set_index(AScolsSubset),
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
            "nucleotides different between the pair. Transcript pairs with 5'/3' "
            "transcript length variation and/or nonoverlapping coordinates "
            "(and no IR or altnerative exon/donor/acceptor) are "
            "counted in the first column with no black dots.".format(
                len(pairAS),
                pairAS["gene_id"].nunique()
            )
        )
    else:
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
                color=sns.color_palette("colorblind", 15).as_hex()[
                    boxCols.index(col)],
            )
    upset.plot()
    plt.subplots_adjust(right=1.00001)
    plt.suptitle(title)


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

# Import trand plotting functions
#from trand import plot_functions as PF


distanceDf = pd.read_csv(
    "C:/Users/knife/OneDrive - University of Florida/McIntyre Lab/MyTranD/mysrc/my_test_csv.csv")

for ref in ["mel", "sim"]:

    print("\n{} Reference:".format(ref))
    # Get union reference UJC flag file
    unionUJCDf = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_UJC_flags.csv".format(ref, ref), low_memory=False)
    unionGeneDf = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_gene_flags.csv".format(ref, ref), low_memory=False)

    # Merge union ref UJC and gene flags
    unionDf = pd.merge(
        unionUJCDf,
        unionGeneDf.drop(columns=["transcript_id"]),
        how="outer",
        on=["gene_id"],
        suffixes=["_in_UJC", "_in_gene"],
        indicator="merge_check",
        validate="m:1"
    )
    if len(unionDf) != unionDf["merge_check"].value_counts()["both"]:
        print("ERROR: Missing some union reference gene_id in input files")

    # Get consolidation key file of reads + union
    combConsolDf = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/consolidate_baseline_w_mapped_union_ref/readsMapRef2{}_transcript_id_2_consolidation_id.csv".format(ref), low_memory=False)

    # Count and flag mel and sim reads consolidated with each union UJC
    combConsolDf["num_mel_read_in_UJC"] = combConsolDf.groupby("consolidation_transcript_id")["transcript_id"].transform(
        lambda x: len(x[x.str.startswith("mel_read_")]))
    combConsolDf["num_sim_read_in_UJC"] = combConsolDf.groupby("consolidation_transcript_id")["transcript_id"].transform(
        lambda x: len(x[x.str.startswith("sim_read_")]))
    combConsolDf["num_union_in_UJC"] = combConsolDf.groupby("consolidation_transcript_id")["transcript_id"].transform(
        lambda x: len(x[x.str.startswith("tr_")]))
    combConsolDf["flag_mel_read_in_UJC"] = np.where(
        combConsolDf["num_mel_read_in_UJC"] > 0,
        1,
        0
    )
    combConsolDf["flag_sim_read_in_UJC"] = np.where(
        combConsolDf["num_sim_read_in_UJC"] > 0,
        1,
        0
    )

    # Flag rows of union reference transcripts
    combConsolDf["flag_in_union"] = np.where(
        combConsolDf["num_union_in_UJC"] > 0,
        1,
        0
    )

    # # Merge in symbols from FB617 annotation file
    # combConsolDf = pd.merge(
    #     combConsolDf,
    #     annotDf[["symbol", "primary_FBgn"]],
    #     how="outer",
    #     left_on="gene_id",
    #     right_on="primary_FBgn",
    #     indicator="merge_check"
    # )

    # keep only rows of union reference transcripts
    unionConsolDf1 = combConsolDf[combConsolDf["transcript_id"].str.startswith(
        "tr_")].copy()
    unionConsolDf = combConsolDf[combConsolDf["transcript_id"].str.startswith(
        "tr_")].copy().drop(columns=["consolidation_transcript_id"])

    # Merge with union reference UJC flag file
    unionMergeDf = pd.merge(
        unionDf.drop(columns=["merge_check"]),
        unionConsolDf,
        how="outer",
        left_on="consolidation_transcript_id",
        suffixes=["_unionRef", "_withReads"],
        right_on="transcript_id",
        indicator="merge_check",
        validate="1:1"
    )
    if len(unionMergeDf) != unionMergeDf["merge_check"].value_counts()["both"]:
        print("ERROR: Missing some union reference UJC transcript_id values")

    # Count number of species-specific UJC in reference have a long read in the opposite species
    print("Of the {0} UJC in both species, {1} ({2:.2%}) are found in the {3} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]),
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] +
            unionMergeDf["flag_sim_read_in_UJC"] == 2]),
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"]+unionMergeDf["flag_sim_read_in_UJC"]
            == 2])/len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]),
        "sim"
    ))
    print("Of the {0} UJC in both species, {1} ({2:.2%}) are found in the {3} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]),
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] +
            unionMergeDf["flag_mel_read_in_UJC"] == 2]),
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"]+unionMergeDf["flag_mel_read_in_UJC"]
            == 2])/len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]),
        "mel"
    ))
    print("Of the {0} UJC in both species, {1} ({2:.2%}) are found in both the {3} AND {4} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]),
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] +
            unionMergeDf["flag_sim_read_in_UJC"]+unionMergeDf["flag_mel_read_in_UJC"] == 3]),
        len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"]+unionMergeDf["flag_sim_read_in_UJC"] +
            unionMergeDf["flag_mel_read_in_UJC"] == 3])/len(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]),
        "mel",
        "sim"
    ))
    print("Of the {0} UJC in {1} only (in shared gene loci), {2} ({3:.2%}) are found in the {4} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "mel",
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_sim_read_in_UJC"] == 3]),
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_sim_read_in_UJC"]
            == 3])/len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "sim"
    ))
    print("Of the {0} UJC in {1} only (in shared gene loci), {2} ({3:.2%}) are found in the {4} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "mel",
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_mel_read_in_UJC"] == 3]),
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_mel_read_in_UJC"]
            == 3])/len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "mel"
    ))
    print("Of the {0} UJC in {1} only (in shared gene loci), {2} ({3:.2%}) are found in the {4} long reads and NOT the {5} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "mel",
        len(unionMergeDf[(unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] +
            unionMergeDf["flag_mel_read_in_UJC"] == 3) & (unionMergeDf["flag_sim_read_in_UJC"] == 0)]),
        len(unionMergeDf[(unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_mel_read_in_UJC"] == 3) & (
            unionMergeDf["flag_sim_read_in_UJC"] == 0)])/len(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "mel",
        "sim"
    ))
    print("Of the {0} UJC in {1} only (in shared gene loci), {2} ({3:.2%}) are found in the {4} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "sim",
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_mel_read_in_UJC"] == 3]),
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_mel_read_in_UJC"]
            == 3])/len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "mel"
    ))
    print("Of the {0} UJC in {1} only (in shared gene loci), {2} ({3:.2%}) are found in the {4} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "sim",
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_sim_read_in_UJC"] == 3]),
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_sim_read_in_UJC"]
            == 3])/len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "sim"
    ))
    print("Of the {0} UJC in {1} only (in shared gene loci), {2} ({3:.2%}) are found in the {4} long reads and NOT the {5} long reads".format(
        len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] +
            unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "sim",
        len(unionMergeDf[(unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] +
            unionMergeDf["flag_sim_read_in_UJC"] == 3) & (unionMergeDf["flag_mel_read_in_UJC"] == 0)]),
        len(unionMergeDf[(unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_sim_read_in_UJC"] == 3) & (
            unionMergeDf["flag_mel_read_in_UJC"] == 0)])/len(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]),
        "sim",
        "mel"
    ))

# Look at shared novel UJC in novel loci of reads
    combConsolGene = combConsolDf.groupby("gene_id")[
        ["flag_in_union", "flag_mel_read_in_UJC", "flag_sim_read_in_UJC"]].max()
    novelGeneDf = combConsolGene[combConsolGene["flag_in_union"] == 0].reset_index(
    )

    # Drop FBgn (!!! SOME FBGN ARE IN HERE BECAUSE SOMEHOW THE GFFREAD GTF TO FA CONVERSINO MISSED THEM)
    novelGeneDf = novelGeneDf[~novelGeneDf["gene_id"].str.startswith("FBgn")]
    novelDf = combConsolDf[combConsolDf["gene_id"].isin(
        novelGeneDf["gene_id"])].copy()

    # Count novel gene loci and UJC and how many are in the reads of mel/sim/both
    print("{} gene loci not in union reference".format(len(novelGeneDf)))
    print("{} novel gene loci in {} reads only, {} novel gene loci in {} reads only, {} novel gene loci in both".format(
        len(novelGeneDf[(novelGeneDf["flag_mel_read_in_UJC"] == 1) & (
            novelGeneDf["flag_sim_read_in_UJC"] == 0)]),
        "mel",
        len(novelGeneDf[(novelGeneDf["flag_sim_read_in_UJC"] == 1) & (
            novelGeneDf["flag_mel_read_in_UJC"] == 0)]),
        "sim",
        len(novelGeneDf[novelGeneDf["flag_mel_read_in_UJC"] +
            novelGeneDf["flag_sim_read_in_UJC"] == 2])
    ))
    print("{} UJC in {} reads, {} UJC in {} reads, {} UJC in either, {} UJC in both".format(
        len(novelDf[novelDf["flag_mel_read_in_UJC"] == 1]),
        "mel",
        len(novelDf[novelDf["flag_sim_read_in_UJC"] == 1]),
        "sim",
        len(novelDf[novelDf["flag_mel_read_in_UJC"] +
            novelDf["flag_sim_read_in_UJC"] > 0]),
        len(novelDf[novelDf["flag_mel_read_in_UJC"] +
            novelDf["flag_sim_read_in_UJC"] == 2])
    ))

    # Get TranD mel vs sim distance output
    melVsimDf = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_mel_vs_sim_2_{}/mel_vs_sim_2_{}_pairwise_transcript_distance.csv".format(ref, ref), low_memory=False)

    # Get CSV of sampleID counts and flags for baseline UJC
    melReadDf = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/baseline_selection_consolidation/mel_2_{}_consolidation_sampleID.csv".format(ref), low_memory=False)
    simReadDf = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/baseline_selection_consolidation/sim_2_{}_consolidation_sampleID.csv".format(ref), low_memory=False)

    # Flag male-limited (in at least one male sample and no female samples)
    melMaleCols = [c for c in melReadDf.columns if "flag" in c and "_m_" in c]
    simMaleCols = [c for c in simReadDf.columns if "flag" in c and "_m_" in c]
    melFemaleCols = [
        c for c in melReadDf.columns if "flag" in c and "_f_" in c]
    simFemaleCols = [
        c for c in simReadDf.columns if "flag" in c and "_f_" in c]

    # Concatenate reads and merge sampleID flags with reads+union
    sampleCombDf = pd.merge(
        pd.concat([melReadDf, simReadDf], ignore_index=True).rename(
            columns={"consolidation_transcript_id": "transcript_id"}).fillna(0),
        combConsolDf,
        how="outer",
        on="transcript_id",
        indicator="merge_check",
        validate="1:1"
    )
    if sampleCombDf["merge_check"].value_counts()["left_only"] != 0:
        print("ERROR: Missing baseline read UJC transcript_id values")

    # Flag male-liimited and female-limited baseline UJC in each species
    sampleCombDf["flag_m_limited_mel"] = np.where(
        (sampleCombDf[melMaleCols].sum(axis=1) > 0)
        & (sampleCombDf[melFemaleCols].sum(axis=1) == 0),
        1,
        0
    )
    sampleCombDf["flag_f_limited_mel"] = np.where(
        (sampleCombDf[melFemaleCols].sum(axis=1) > 0)
        & (sampleCombDf[melMaleCols].sum(axis=1) == 0),
        1,
        0
    )
    sampleCombDf["flag_both_sexes_mel"] = np.where(
        (sampleCombDf[melFemaleCols].sum(axis=1) > 0)
        & (sampleCombDf[melMaleCols].sum(axis=1) > 0),
        1,
        0
    )
    sampleCombDf["flag_m_limited_sim"] = np.where(
        (sampleCombDf[simMaleCols].sum(axis=1) > 1)
        & (sampleCombDf[simFemaleCols].sum(axis=1) == 0),
        1,
        0
    )
    sampleCombDf["flag_f_limited_sim"] = np.where(
        (sampleCombDf[simFemaleCols].sum(axis=1) > 1)
        & (sampleCombDf[simMaleCols].sum(axis=1) == 0),
        1,
        0
    )
    sampleCombDf["flag_both_sexes_sim"] = np.where(
        (sampleCombDf[simFemaleCols].sum(axis=1) > 0)
        & (sampleCombDf[simMaleCols].sum(axis=1) > 0),
        1,
        0
    )
    sampleCols = ['flag_mel_read_in_UJC', 'flag_sim_read_in_UJC', 'flag_in_union',
                  'flag_m_limited_mel', 'flag_f_limited_mel', 'flag_both_sexes_mel',
                  'flag_m_limited_sim', 'flag_f_limited_sim', 'flag_both_sexes_sim']
    consolSampleComb = sampleCombDf.groupby(["gene_id", "consolidation_transcript_id"])[
        sampleCols].max().reset_index()
    consolSampleCombWunion = pd.merge(
        consolSampleComb,
        unionConsolDf1[["consolidation_transcript_id", "transcript_id"]],
        how="outer",
        on="consolidation_transcript_id",
        indicator="merge_check",
        validate="1:m"
    )
    unionSample = consolSampleCombWunion[consolSampleCombWunion["merge_check"] == "both"].copy(
    )
    geneSampleComb = consolSampleComb.groupby(
        "gene_id")[sampleCols].max().reset_index()

    # Count sex-limited UJC in the reads
    print("{0} of the {1} annotated UJC in the {2} long reads ({3:.2%}) are male-limited".format(
        consolSampleComb[consolSampleComb["flag_in_union"]
                         == 1]["flag_m_limited_mel"].sum(),
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_mel_read_in_UJC"] == 2]),
        "mel",
        consolSampleComb[consolSampleComb["flag_in_union"] == 1]["flag_m_limited_mel"].sum(
        )/len(consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_mel_read_in_UJC"] == 2]),
    ))
    print("{0} of the {1} annotated UJC in the {2} long reads ({3:.2%}) are female-limited".format(
        consolSampleComb[consolSampleComb["flag_in_union"]
                         == 1]["flag_f_limited_mel"].sum(),
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_mel_read_in_UJC"] == 2]),
        "mel",
        consolSampleComb[consolSampleComb["flag_in_union"] == 1]["flag_f_limited_mel"].sum(
        )/len(consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_mel_read_in_UJC"] == 2]),
    ))
    print("{0} of the {1} genes with annotated UJC in the {2} long reads ({3:.2%}) have at least one sex-limited isoform".format(
        consolSampleComb[(consolSampleComb["flag_in_union"] == 1) & (
            consolSampleComb["flag_m_limited_mel"]+consolSampleComb["flag_f_limited_mel"] > 0)]["gene_id"].nunique(),
        consolSampleComb[consolSampleComb["flag_in_union"] +
                         consolSampleComb["flag_mel_read_in_UJC"] == 2]["gene_id"].nunique(),
        "mel",
        consolSampleComb[(consolSampleComb["flag_in_union"] == 1) & (consolSampleComb["flag_m_limited_mel"]+consolSampleComb["flag_f_limited_mel"] > 0)
                         ]["gene_id"].nunique()/consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_mel_read_in_UJC"] == 2]["gene_id"].nunique(),
    ))
    print("{0} of the {1} annotated UJC in the {2} long reads ({3:.2%}) are male-limited".format(
        consolSampleComb[consolSampleComb["flag_in_union"]
                         == 1]["flag_m_limited_sim"].sum(),
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_sim_read_in_UJC"] == 2]),
        "sim",
        consolSampleComb[consolSampleComb["flag_in_union"] == 1]["flag_m_limited_sim"].sum(
        )/len(consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_sim_read_in_UJC"] == 2]),
    ))
    print("{0} of the {1} annotated UJC in the {2} long reads ({3:.2%}) are female-limited".format(
        consolSampleComb[consolSampleComb["flag_in_union"]
                         == 1]["flag_f_limited_sim"].sum(),
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_sim_read_in_UJC"] == 2]),
        "sim",
        consolSampleComb[consolSampleComb["flag_in_union"] == 1]["flag_f_limited_sim"].sum(
        )/len(consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_sim_read_in_UJC"] == 2]),
    ))
    print("{0} of the {1} genes with annotated UJC in the {2} long reads ({3:.2%}) have at least one sex-limited isoform".format(
        consolSampleComb[(consolSampleComb["flag_in_union"] == 1) & (
            consolSampleComb["flag_m_limited_sim"]+consolSampleComb["flag_f_limited_sim"] > 0)]["gene_id"].nunique(),
        consolSampleComb[consolSampleComb["flag_in_union"] +
                         consolSampleComb["flag_sim_read_in_UJC"] == 2]["gene_id"].nunique(),
        "sim",
        consolSampleComb[(consolSampleComb["flag_in_union"] == 1) & (consolSampleComb["flag_m_limited_sim"]+consolSampleComb["flag_f_limited_sim"] > 0)
                         ]["gene_id"].nunique()/consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_sim_read_in_UJC"] == 2]["gene_id"].nunique(),
    ))
    print("{0} of the {1} annotated UJC in the both species long reads ({2:.2%}) are male-limited in both species".format(
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_m_limited_mel"]+consolSampleComb["flag_m_limited_sim"] == 3]),
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_mel_read_in_UJC"]+consolSampleComb["flag_sim_read_in_UJC"] == 3]),
        len(consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_m_limited_mel"]+consolSampleComb["flag_m_limited_sim"] == 3])/len(
            consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_mel_read_in_UJC"]+consolSampleComb["flag_sim_read_in_UJC"] == 3]),
    ))
    print("{0} of the {1} annotated UJC in the both species long reads ({2:.2%}) are female-limited in both species".format(
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_f_limited_mel"]+consolSampleComb["flag_f_limited_sim"] == 3]),
        len(consolSampleComb[consolSampleComb["flag_in_union"] +
            consolSampleComb["flag_mel_read_in_UJC"]+consolSampleComb["flag_sim_read_in_UJC"] == 3]),
        len(consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_f_limited_mel"]+consolSampleComb["flag_f_limited_sim"] == 3])/len(
            consolSampleComb[consolSampleComb["flag_in_union"]+consolSampleComb["flag_mel_read_in_UJC"]+consolSampleComb["flag_sim_read_in_UJC"] == 3]),
    ))

    # Output transcript pairs for each species of transcripts that are:
    #   1) male-limited and female-limited
    #   2) male-limited and both
    #   3) female-limited and both
    # Get union reference pairwise distance for sex limited subsets
    unionPairDf1 = pd.read_csv(
        "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_pairwise_transcript_distance.csv".format(ref, ref), low_memory=False)

    if ref == "mel":
        annotDf = pd.read_csv(
            "/Volumes/blue/mcintyre/share/references/dmel_fb617/fbgn_annotation_ID.csv", low_memory=False)
    # Merge in symbols from FB617 annotation file
    unionPairDf = pd.merge(
        unionPairDf1,
        annotDf[["symbol", "primary_FBgn"]],
        how="outer",
        left_on="gene_id",
        right_on="primary_FBgn",
        indicator="merge_check"
    )

    subset1 = unionPairDf[
        ((unionPairDf["transcript_1"].isin(unionSample[unionSample["flag_m_limited_"+ref] == 1]["transcript_id"])) &
         (unionPairDf["transcript_2"].isin(unionSample[unionSample["flag_f_limited_"+ref] == 1]["transcript_id"])))
        | ((unionPairDf["transcript_1"].isin(unionSample[unionSample["flag_f_limited_"+ref] == 1]["transcript_id"])) & (unionPairDf["transcript_2"].isin(unionSample[unionSample["flag_m_limited_"+ref] == 1]["transcript_id"])))
    ]
    subset1.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_m_f_limited_{}_pairwise_transcript_distance.csv".format(ref, ref, ref), index=False)
    subset2 = unionPairDf[
        ((unionPairDf["transcript_1"].isin(unionSample[unionSample["flag_m_limited_"+ref] == 1]["transcript_id"])) &
         (unionPairDf["transcript_2"].isin(unionSample[unionSample["flag_both_sexes_"+ref] == 1]["transcript_id"])))
        | ((unionPairDf["transcript_1"].isin(unionSample[unionSample["flag_both_sexes_"+ref] == 1]["transcript_id"])) & (unionPairDf["transcript_2"].isin(unionSample[unionSample["flag_m_limited_"+ref] == 1]["transcript_id"])))
    ]
    subset2.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_m_limited_bothSex_{}_pairwise_transcript_distance.csv".format(ref, ref, ref), index=False)
    subset3 = unionPairDf[
        ((unionPairDf["transcript_1"].isin(unionSample[unionSample["flag_both_sexes_"+ref] == 1]["transcript_id"]))
         & (unionPairDf["transcript_2"].isin(unionSample[unionSample["flag_f_limited_"+ref] == 1]["transcript_id"])))
        | ((unionPairDf["transcript_1"].isin(unionSample[unionSample["flag_f_limited_"+ref] == 1]["transcript_id"])) & (unionPairDf["transcript_2"].isin(unionSample[unionSample["flag_both_sexes_"+ref] == 1]["transcript_id"])))
    ]
    subset2.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_f_limited_bothSex_{}_pairwise_transcript_distance.csv".format(ref, ref, ref), index=False)

    # Plot all pairs for each set
#    outdir = "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}".format(ref)
    outdir = "/Users/adalena/mclab/SHARE/McIntyre_Lab/Transcript_orthologs/figures/union2{}_sex_specific_upset".format(
        ref)
    plot_pair_AS_upset_nt_box(subset1, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir, "union2{}_m_f_limited_{}".format(ref, ref)), drop_5_3_no_shared_nt=True)
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir,
                "union2{}_m_f_limited".format(ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subset2, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir, "union2{}_m_limited_bothSex_{}".format(ref, ref)), drop_5_3_no_shared_nt=True)
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir,
                "union2{}_m_limited_bothSex_{}".format(ref, ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subset3, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir, "union2{}_f_limited_bothSex_{}".format(ref, ref)), drop_5_3_no_shared_nt=True)
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir,
                "union2{}_f_limited_bothSex_{}".format(ref, ref)), dpi=600, format="png")

    # Make the above sex-specific reference plots for the following genes
    if ref == "mel":
        # Sxl (FBgn0264270)
        # fru (FBgn0004652)
        # dsx (FBgn0000504)
        # Found weird transcripts that need to be dropped
        dropList = ["tr_FBgn0004652_1", "tr_FBgn0004652_2"]
        for symbol, gene in [("Sxl", "FBgn0264270"), ("fru", "FBgn0004652"), ("dsx", "FBgn0000504")]:
            geneSubset1 = subset1[(subset1["gene_id"] == gene) & (
                ~subset1["transcript_1"].isin(dropList)) & (~subset1["transcript_2"].isin(dropList))]
            geneSubset2 = subset2[(subset2["gene_id"] == gene) & (
                ~subset2["transcript_1"].isin(dropList)) & (~subset2["transcript_2"].isin(dropList))]
            geneSubset3 = subset3[(subset3["gene_id"] == gene) & (
                ~subset3["transcript_1"].isin(dropList)) & (~subset3["transcript_2"].isin(dropList))]
            if len(geneSubset1) > 0:
                plot_pair_AS_upset_nt_box(geneSubset1, "{}/{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.rtf".format(
                    outdir, symbol, "union2{}_m_f_limited_{}".format(ref, ref), symbol))
                plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.png".format(
                    outdir, "union2{}_m_f_limited_".format(ref), symbol), dpi=600, format="png")
            if len(geneSubset2) > 0:
                plot_pair_AS_upset_nt_box(geneSubset2, "{}/{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.rtf".format(
                    outdir, symbol, "union2{}_m_limited_bothSex_{}".format(ref, ref), symbol))
                plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.png".format(outdir,
                            "union2{}_m_limited_bothSex_{}".format(ref, ref), symbol), dpi=600, format="png")
            if len(geneSubset3) > 0:
                plot_pair_AS_upset_nt_box(geneSubset3, "{}/{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.rtf".format(
                    outdir, symbol, "union2{}_f_limited_bothSex_{}".format(ref, ref), symbol))
                plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.png".format(outdir,
                            "union2{}_f_limited_bothSex_{}".format(ref, ref), symbol), dpi=600, format="png")

    # Output transcript pairs from union reference of transcripts that are:
    #   1) mel-specific and sim-specific
    #   2) mel-specific and both
    #   3) sim-specific and both
    # Get union reference pairwise distance for species-specific subsets
    subsetMelSim1 = unionPairDf[
        ((unionPairDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] == 1]["consolidation_transcript_id"])) & (
            unionPairDf["transcript_2"].isin(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] == 1]["consolidation_transcript_id"])))
        | ((unionPairDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] == 1]["consolidation_transcript_id"])) & (unionPairDf["transcript_2"].isin(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] == 1]["consolidation_transcript_id"])))
    ]
    subsetMelSim1.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_mel_sim_specific_{}_pairwise_transcript_distance.csv".format(ref, ref, ref), index=False)
    subsetMelSim2 = unionPairDf[
        ((unionPairDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] == 1]["consolidation_transcript_id"])) & (
            unionPairDf["transcript_2"].isin(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]["consolidation_transcript_id"])))
        | ((unionPairDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]["consolidation_transcript_id"])) & (unionPairDf["transcript_2"].isin(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] == 1]["consolidation_transcript_id"])))
    ]
    subsetMelSim2.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_mel_specific_bothSpp_{}_pairwise_transcript_distance.csv".format(ref, ref, ref), index=False)
    subsetMelSim3 = unionPairDf[
        ((unionPairDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]["consolidation_transcript_id"])) & (
            unionPairDf["transcript_2"].isin(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] == 1]["consolidation_transcript_id"])))
        | ((unionPairDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] == 1]["consolidation_transcript_id"])) & (unionPairDf["transcript_2"].isin(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"] == 1]["consolidation_transcript_id"])))
    ]
    subsetMelSim3.to_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}/union2{}_sim_specific_bothSpp_{}_pairwise_transcript_distance.csv".format(ref, ref, ref), index=False)

    # Plot all pairs for each set
    outdir = "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/reference_comparison/TranD_consol_union_2_{}".format(
        ref)
    plot_pair_AS_upset_nt_box(subsetMelSim1, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir, "union2{}_mel_sim_specific_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir,
                "union2{}_mel_sim_specific_".format(ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subsetMelSim2, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir, "union2{}_mel_specific_bothSpp_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir,
                "union2{}_mel_specific_bothSpp_{}".format(ref, ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subsetMelSim3, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir, "union2{}_sim_specific_bothSpp_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir,
                "union2{}_sim_specific_bothSpp_{}".format(ref, ref)), dpi=600, format="png")

    # Make the above union reference plots for the following genes
    if ref == "mel":
        # Sxl (FBgn0264270)
        # fru (FBgn0004652)
        # dsx (FBgn0000504)
        # Found weird transcripts that need to be dropped
        dropList = ["tr_FBgn0004652_1", "tr_FBgn0004652_2"]
        for symbol, gene in [("Sxl", "FBgn0264270"), ("fru", "FBgn0004652"), ("dsx", "FBgn0000504")]:
            geneMelSim1 = subsetMelSim1[(subsetMelSim1["gene_id"] == gene) & (
                ~subsetMelSim1["transcript_1"].isin(dropList)) & (~subsetMelSim1["transcript_2"].isin(dropList))]
            geneMelSim2 = subsetMelSim2[(subsetMelSim2["gene_id"] == gene) & (
                ~subsetMelSim2["transcript_1"].isin(dropList)) & (~subsetMelSim2["transcript_2"].isin(dropList))]
            geneMelSim3 = subsetMelSim3[(subsetMelSim3["gene_id"] == gene) & (
                ~subsetMelSim3["transcript_1"].isin(dropList)) & (~subsetMelSim3["transcript_2"].isin(dropList))]
            if len(geneMelSim1) > 0:
                plot_pair_AS_upset_nt_box(geneMelSim1, "{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.rtf".format(
                    outdir, "union2{}_mel_sim_specific_{}".format(ref, ref), symbol))
                plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.png".format(
                    outdir, "union2{}_mel_sim_specific_".format(ref), symbol), dpi=600, format="png")
            if len(geneMelSim2) > 0:
                plot_pair_AS_upset_nt_box(geneMelSim2, "{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.rtf".format(
                    outdir, "union2{}_mel_specific_bothSpp_{}".format(ref, ref), symbol))
                plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.png".format(outdir,
                            "union2{}_mel_specific_bothSpp_{}".format(ref, ref), symbol), dpi=600, format="png")
            if len(geneMelSim3) > 0:
                plot_pair_AS_upset_nt_box(geneMelSim3, "{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.rtf".format(
                    outdir, "union2{}_sim_specific_bothSpp_{}".format(ref, ref), symbol))
                plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot_{}.png".format(outdir,
                            "union2{}_sim_specific_bothSpp_{}".format(ref, ref), symbol), dpi=600, format="png")

    # Get TranD pairwise distance CSV of reads vs. union ref
    distDf = pd.read_csv("/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/consolidate_baseline_w_mapped_union_ref/TranD_{}_union_vs_reads/union_vs_reads_2_{}_pairwise_transcript_distance.csv".format(ref, ref), low_memory=False)
    distDf["transcript_1"] = distDf["transcript_1"].str[:-6]
    distDf["transcript_2"] = distDf["transcript_2"].str[:-5]

    # Subset for species-specific UJC vs. reads in that species
    subsetReadMelSim1 = distDf[
        (distDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_mel_only_in_UJC"] +
         unionMergeDf["flag_in_both_species_in_gene"] == 2]["consolidation_transcript_id"]))
        & (distDf["flag_min_match_union"] == 1)]
    subsetReadMelSim2 = distDf[
        (distDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_sim_only_in_UJC"] +
         unionMergeDf["flag_in_both_species_in_gene"] == 2]["consolidation_transcript_id"]))
        & (distDf["flag_min_match_union"] == 1)]
    subsetReadMelSim3 = distDf[(
        distDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] == 2]["consolidation_transcript_id"]))
        & (distDf["flag_min_match_union"] == 1)]

    # Plot all pairs for each set
    outdir2 = "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/consolidate_baseline_w_mapped_union_ref/TranD_{}_union_vs_reads".format(
        ref)
    plot_pair_AS_upset_nt_box(subsetReadMelSim1, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir2, "union2{}_reads_mel_specific_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir2,
                "union2{}_mel_sim_specific_".format(ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subsetReadMelSim2, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir2, "union2{}_reads_sim_specific_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir2,
                "union2{}_mel_specific_bothSpp_{}".format(ref, ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subsetReadMelSim3, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir2, "union2{}_reads_bothSpp_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir2,
                "union2{}_sim_specific_bothSpp_{}".format(ref, ref)), dpi=600, format="png")

    # Subset for species-specific UJC vs. reads for UJC with reads in the species-specific and not in the other species
    subsetReadMelSim4 = distDf[
        (distDf["transcript_1"].isin(unionMergeDf[(unionMergeDf["flag_in_mel_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] +
         unionMergeDf["flag_mel_read_in_UJC"] == 3) & (unionMergeDf["flag_sim_read_in_UJC"] == 0)]["consolidation_transcript_id"]))
        & (distDf["flag_min_match_union"] == 1)]
    subsetReadMelSim5 = distDf[
        (distDf["transcript_1"].isin(unionMergeDf[(unionMergeDf["flag_in_sim_only_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"] +
         unionMergeDf["flag_sim_read_in_UJC"] == 3) & (unionMergeDf["flag_mel_read_in_UJC"] == 0)]["consolidation_transcript_id"]))
        & (distDf["flag_min_match_union"] == 1)]
    subsetReadMelSim6 = distDf[(
        distDf["transcript_1"].isin(unionMergeDf[unionMergeDf["flag_in_both_species_in_UJC"]+unionMergeDf["flag_in_both_species_in_gene"]+unionMergeDf["flag_mel_read_in_UJC"]+unionMergeDf["flag_sim_read_in_UJC"] == 4]["consolidation_transcript_id"]))
        & (distDf["flag_min_match_union"] == 1)]

    # Plot all pairs for each set
    outdir2 = "/Volumes/blue/mcintyre/share/transcript_distance/dros_analysis/consolidate_baseline_w_mapped_union_ref/TranD_{}_union_vs_reads".format(
        ref)
    plot_pair_AS_upset_nt_box(subsetReadMelSim4, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir2, "union2{}_readInMelnotSim_mel_specific_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir2,
                "union2{}_readInMelnotSim_mel_specific_".format(ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subsetReadMelSim5, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir2, "union2{}_readInSimnotMel_sim_specific_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir2,
                "union2{}_readInSimnotMel_sim_specific_{}".format(ref, ref)), dpi=600, format="png")
    plot_pair_AS_upset_nt_box(subsetReadMelSim6, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(
        outdir2, "union2{}_readInBoth_bothSpp_{}".format(ref, ref)))
    plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir2,
                "union2{}_readInBoth_bothSpp_{}".format(ref, ref)), dpi=600, format="png")


# OUTPUT OF SCRIPT:
    # mel Reference:
    # Of the 15948 UJC in both species, 10935 (68.57%) are found in the sim long reads
    # Of the 15948 UJC in both species, 10663 (66.86%) are found in the mel long reads
    # Of the 15948 UJC in both species, 9654 (60.53%) are found in both the mel AND sim long reads
    # Of the 11040 UJC in mel only (in shared gene loci), 2686 (24.33%) are found in the sim long reads
    # Of the 11040 UJC in mel only (in shared gene loci), 3334 (30.20%) are found in the mel long reads
    # Of the 11040 UJC in mel only (in shared gene loci), 1295 (11.73%) are found in the mel long reads and NOT the sim long reads
    # Of the 8411 UJC in sim only (in shared gene loci), 892 (10.61%) are found in the mel long reads
    # Of the 8411 UJC in sim only (in shared gene loci), 1524 (18.12%) are found in the sim long reads
    # Of the 8411 UJC in sim only (in shared gene loci), 868 (10.32%) are found in the sim long reads and NOT the mel long reads

    # sim Reference:
    # Of the 15372 UJC in both species, 10581 (68.83%) are found in the sim long reads
    # Of the 15372 UJC in both species, 10270 (66.81%) are found in the mel long reads
    # Of the 15372 UJC in both species, 9323 (60.65%) are found in both the mel AND sim long reads
    # Of the 12014 UJC in mel only (in shared gene loci), 1850 (15.40%) are found in the sim long reads
    # Of the 12014 UJC in mel only (in shared gene loci), 2567 (21.37%) are found in the mel long reads
    # Of the 12014 UJC in mel only (in shared gene loci), 1147 (9.55%) are found in the mel long reads and NOT the sim long reads
    # Of the 8747 UJC in sim only (in shared gene loci), 1429 (16.34%) are found in the mel long reads
    # Of the 8747 UJC in sim only (in shared gene loci), 2126 (24.31%) are found in the sim long reads
    # Of the 8747 UJC in sim only (in shared gene loci), 1019 (11.65%) are found in the sim long reads and NOT the mel long reads

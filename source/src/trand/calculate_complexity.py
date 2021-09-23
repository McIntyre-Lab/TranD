#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 16:46:40 2021

@author: adalena.nanni
"""

import pandas as pd
import matplotlib.pyplot as plt

# Import plot functions
from . import plot_functions as PF


def calculate_complexity(outdir, gtf_df, skip_plots, prefix=None):
    """
    Calculate complexity measures from gtf dataframe:
        1) Transcripts per gene
        2) Unique exons per gene
        3) Exons per transcript
    """
    # Get transcript and gene level dataframes using dataframe from gtf
    gtf_df["num_exon"] = 1
    gtf_df["exon_id"] = (
        gtf_df["seqname"].map(str)
        + ":"
        + gtf_df["start"].map(str)
        + ":"
        + gtf_df["end"].map(str)
        + ":"
        + gtf_df["strand"].map(str)
    )
    transcriptDF = (
        gtf_df.groupby("transcript_id").agg({"num_exon": "sum"}).reset_index()
    )
    geneDF = (
        gtf_df.groupby("gene_id")
        .agg({"transcript_id": "nunique", "num_exon": "sum", "exon_id": "nunique"})
        .reset_index()
        .rename(columns={"transcript_id": "num_transcript", "exon_id": "num_uniq_exon"})
    )

    counts = pd.DataFrame()
    # Get number of transcripts
    counts.loc[0, "num_transcript"] = geneDF["num_transcript"].sum()

    # Get number of genes
    counts.loc[0, "num_gene"] = len(geneDF)

    # Get number of exons (same exon in multiple transcripts counted multiple times)
    counts.loc[0, "num_exon"] = geneDF["num_exon"].sum()

    # Get number of unique exons (can still overlap but not have the same coordinates)
    # Must use GTF and not gene level counts since unique exons can be in multiple genes
    counts.loc[0, "num_uniqExon"] = len(
        gtf_df[["seqname", "start", "end", "strand"]].drop_duplicates()
    )

    # Get minimum values
    counts.loc[0, "min_transcriptPerGene"] = geneDF["num_transcript"].min()
    counts.loc[0, "min_exonPerTranscript"] = transcriptDF["num_exon"].min()
    counts.loc[0, "min_exonPerGene"] = geneDF["num_uniq_exon"].min()

    # Get Q1 values
    counts.loc[0, "q1_transcriptPerGene"] = geneDF["num_transcript"].quantile(0.25)
    counts.loc[0, "q1_exonPerTranscript"] = transcriptDF["num_exon"].quantile(0.25)
    counts.loc[0, "q1_exonPerGene"] = geneDF["num_uniq_exon"].quantile(0.25)

    # Get median values
    counts.loc[0, "median_transcriptPerGene"] = geneDF["num_transcript"].median()
    counts.loc[0, "median_exonPerTranscript"] = transcriptDF["num_exon"].median()
    counts.loc[0, "median_exonPerGene"] = geneDF["num_uniq_exon"].median()

    # Get Q3 values
    counts.loc[0, "q3_transcriptPerGene"] = geneDF["num_transcript"].quantile(0.75)
    counts.loc[0, "q3_exonPerTranscript"] = transcriptDF["num_exon"].quantile(0.75)
    counts.loc[0, "q3_exonPerGene"] = geneDF["num_uniq_exon"].quantile(0.75)

    # Get max values
    counts.loc[0, "max_transcriptPerGene"] = geneDF["num_transcript"].max()
    counts.loc[0, "max_exonPerTranscript"] = transcriptDF["num_exon"].max()
    counts.loc[0, "max_exonPerGene"] = geneDF["num_uniq_exon"].max()

    # Get mean values
    counts.loc[0, "mean_transcriptPerGene"] = geneDF["num_transcript"].mean()
    counts.loc[0, "mean_exonPerTranscript"] = transcriptDF["num_exon"].mean()
    counts.loc[0, "mean_exonPerGene"] = geneDF["num_uniq_exon"].mean()

    # Get std
    counts.loc[0, "std_transcriptPerGene"] = geneDF["num_transcript"].std()
    counts.loc[0, "std_exonPerTranscript"] = transcriptDF["num_exon"].std()
    counts.loc[0, "std_exonPerGene"] = geneDF["num_uniq_exon"].std()

    # Plot boxplots of complexity measures and output counts to file
    if prefix is None:
        if not skip_plots:
            PF.plot_complexity_box(
                geneDF, transcriptDF, outdir, "{}/complexity_plots.rtf".format(outdir)
            )
            plt.savefig("{}/complexity_plots.png".format(outdir), dpi=600, format="png")
            plt.clf()
        counts.to_csv(
            "{}/transcriptome_complexity_counts.csv".format(outdir), index=False
        )
        return geneDF[["gene_id", "num_uniq_exon"]]
    #        plot_complexity_violin(geneDF,transcriptDF,outdir,legendOut)
    #        plt.savefig("{}/complexity_plots_violin.png".format(outdir),dpi=600,format="png")
    else:
        if not skip_plots:
            PF.plot_complexity_box(
                geneDF,
                transcriptDF,
                outdir,
                "{}/{}_complexity_plots.rtf".format(outdir, prefix),
            )
            plt.savefig(
                "{}/{}_complexity_plots.png".format(outdir, prefix),
                dpi=600,
                format="png",
            )
            plt.clf()
        counts.to_csv(
            "{}/{}_transcriptome_complexity_counts.csv".format(outdir, prefix),
            index=False,
        )

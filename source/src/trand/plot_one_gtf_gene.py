#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 16:13:11 2021

@author: adalena.nanni
"""

import matplotlib.pyplot as plt
from loguru import logger

# Import plot functions
from . import plot_functions as PF

# Import transcriptome summary plot functions
from . import plot_transcriptome as PT


def plot_one_gtf_gene(er_data, ef_data, ir_data, uniqex_data, outdir):
    """
    Plot 1 GTF Pairwise output including the following:

        1) Upset plot for the number of genes with each kind of
            alternative splicing (AS) and the distribution of the average
            nt differences between transcript pairs within each gene
    """
# Plot 1
    # Summary of structural elements of transcriptome, grouped by number of
    #   transcripts per gene
    PT.plot_transcriptome(er_data, ef_data, ir_data, uniqex_data, outdir)
    plt.clf()
    # Make legend for plot 1
    # !!! Move this to plot functions when plot_transcriptome is moved
    legendOut = "{}/transcriptome_summary_plot.rtf".format(outdir)
    legendText = (
        "Summary of structural elements from input GTF. Structural elements "
        "include the frequency of genes within each range of transcript per "
        "gene values (top row, histogram). For each range of transcripts per "
        "gene, the distribution of exon regions per gene (second row, "
        "box plot) represent the number of exons disregarding alternative "
        "donors and acceptors, the number of unique exons per gene (third row, "
        "box plot) represents alternative donor/acceptors and alternative "
        "exons, the proportion of transcripts with intron retention "
        "(fourth row, heatmap), the proportion of exon regions with alternate "
        "donors/acceptors (fifth row, heatmap), the proportion of variable "
        "exon regions (sixth row, heatmap) represents the alternative exons "
        "disregarding alternative donors and acceptors, and the proportion of "
        "variable nucleotides (bottom row, heatmap) represents nucleotides "
        "that are contained in alternative exons and alternative donors and "
        "acceptors."
    )
    with open(legendOut, "w") as outFile:
        PF.start_rtf(outFile)
        outFile.write(
            r"\b Figure. Transcriptome summary \b0 \line {}".format(
                legendText
            )
        )
        PF.end_rtf(outFile)

# Plot 2
    # Distribution of proportion of nucleotide variability across all genes
    try:
        PF.plot_gene_prop_nt_variablility(ef_data,"{}/all_gene_prop_nt_variablility.rtf".format(outdir))
        plt.savefig("{}/all_gene_prop_nt_variablility.png".format(outdir), dpi=600, format="png")
        plt.clf()
    except ValueError as e:
        logger.error("Not enough data to plot nt variability {}".format(e))
        return

# Plot 3
    # Distribution of proportion of nucleotide variability across multi-transcript
    #   genes
    PF.plot_gene_prop_nt_variablility(ef_data,
                                      "{}/multi_xcrpt_gene_prop_nt_variablility.rtf".format(outdir),
                                      multitranscript=True)
    plt.savefig("{}/multi_xcrpt_gene_prop_nt_variablility.png".format(outdir), dpi=600, format="png")
    plt.clf()

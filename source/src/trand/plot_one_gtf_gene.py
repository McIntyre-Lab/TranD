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

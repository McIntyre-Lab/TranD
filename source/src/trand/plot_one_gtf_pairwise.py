#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 16:13:11 2021

@author: adalena.nanni
"""

import matplotlib.pyplot as plt

# Import plot functions
from . import plot_functions as PF


def plot_one_gtf_pairwise(outdir, td_data, prefix):
    """
    Plot 1 GTF Pairwise output including the following:

        1) Upset plot for the number of genes with each kind of alternative splicing (AS)

        2) Upset plot the number of transcript pairs with each kind of AS
           for a single transcriptome with added box plots
           of the number of nt different and the proportion of nt different
    """
# Plot 1
    # Upset plot for the number of genes with each kind of alternative splicing (AS)
    if prefix is not None:
        PF.plot_gene_AS_upset(td_data, "{}/{}_gene_AS_upset.rtf".format(outdir, prefix))
        plt.savefig("{}/{}_gene_AS_upset.png".format(outdir, prefix), dpi=600, format="png")
    else:
        PF.plot_gene_AS_upset(td_data, "{}/gene_AS_upset.rtf".format(outdir))
        plt.savefig("{}/gene_AS_upset.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 2
    # Upset plot the number of transcript pairs with each kind of
    #   AS for a single transcriptome with added box plots
    #   of the number of nt different and the proportion of nt different
    if prefix is not None:
        PF.plot_pair_AS_upset_nt_box(td_data, "{}/{}_transcript_pair_AS_upset_nt_boxplot.rtf".format(outdir, prefix))
        plt.savefig("{}/{}_transcript_pair_AS_upset_nt_boxplot.png".format(outdir, prefix), dpi=600, format="png")
    else:
        PF.plot_pair_AS_upset_nt_box(td_data, "{}/transcript_pair_AS_upset_nt_boxplot.rtf".format(outdir))
        plt.savefig("{}/transcript_pair_AS_upset_nt_boxplot.png".format(outdir), dpi=600, format="png")
    plt.clf()
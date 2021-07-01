#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 16:13:11 2021

@author: adalena.nanni
"""

import matplotlib.pyplot as plt

# Import plot functions
from . import plot_functions as PF


def plot_one_gtf_pairwise(outdir, td_data):
    """
    Plot 1 GTF Pairwise output including the following:

        1) Upset plot for the number of genes with each kind of
            alternative splicing (AS) and the distribution of the average
            nt differences between transcript pairs within each gene
    """
# ##### Plot 1 ######
    # Upset plot for the number of genes with each kind of alternative splicing (AS)
    #   and the distribution of the average nt differences between transcript pairs
    #   within each gene
    PF.plot_gene_AS_upset_nt_box(td_data, "{}/gene_AS_upset_nt_boxplot.rtf".format(outdir))
    plt.savefig("{}/gene_AS_upset_nt_boxplot.png".format(outdir), dpi=600, format="png")
    plt.clf()

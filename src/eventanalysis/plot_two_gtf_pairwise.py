#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:33:04 2021

@author: adalena.nanni
"""

import matplotlib.pyplot as plt

# Import plot functions
import plot_functions as PF

def plot_two_gtf_pairwise(outdir,md_data,f1_odds,f2_odds,name1,name2):
    # Plot pie chart for number of genes where:
    #   the numbers of transcripts are the same (Equal) between the two datasets,
    #   where one dataset has more transcripts that the other,
    #   or where one dataset has transcripts and the other has none for the gene
    PF.plot_transcript_in_gene_pie(md_data,f1_odds,f2_odds,name1,name2)
#    plt.savefig("test_transcript_in_gene_pie.png",dpi=600,format="png")
    plt.savefig("{}/transcript_in_gene_pie.png".format(outdir),dpi=600,format="png")
    plt.clf()
    
    # Plot pie charts split by 1) all genes, 2) single transcript in at least 1 dataset,
    #   3) multi-transcript in at least 1 dataset, and 4) multi-transcript in both datasets
    PF.plot_transcript_in_gene_split_pie(md_data,f1_odds,f2_odds,name1,name2)
#    plt.savefig("test_transcript_in_gene_split_pie.png",dpi=600,format="png")
    plt.savefig("{}/transcript_in_gene_split_pie.png".format(outdir),dpi=600,format="png")
    plt.clf()
    
    # Upset plot for number of transcripts in genes
    PF.plot_transcript_in_gene_upset(md_data,f1_odds,f2_odds,name1,name2)
#    plt.savefig("test_transcript_in_gene_upset.png",dpi=600,format="png")
    plt.savefig("{}/transcript_in_gene_upset.png".format(outdir),dpi=600,format="png")
    plt.clf()
    
    # Plot stacked bar chart of genes with reciprocal minimum pairs for:
    #   1) all transcripts possible (match or subset)
    #   or 2) only a subset of transcripts possible (partial match or partial subset)
    PF.plot_gene_stack(md_data,name1,name2)
#    plt.savefig("test_transcript_in_gene_stackCount.png",dpi=600,format="png")
    plt.savefig("{}/transcript_in_gene_stackCount.png".format(outdir),dpi=600,format="png")
    plt.clf()
    
    # With percentage of genes instead of count
    PF.plot_gene_stack(md_data,name1,name2,useProp=True)
#    plt.savefig("test_transcript_in_gene_stackProp.png",dpi=600,format="png")
    plt.savefig("{}/transcript_in_gene_stack_proportion.png".format(outdir),dpi=600,format="png")
    plt.clf()    
    
    # For all reciprocal minimum pair transcripts:
    # Upset plot for the number of transcripts with each kind of AS
    PF.plot_recip_min_pair_AS_upset(md_data)
#    plt.savefig("test_recip_min_pair_AS_upset.png",dpi=600,format="png")
    plt.savefig("{}/recip_min_pair_AS_upset.png".format(outdir),dpi=600,format="png")
    plt.clf()

    # For all reciprocal minimum pair transcripts:
    # Upset plot for the number of transcripts with each kind of AS and the distribution
    #   of nt differences between the pairs plotted above
    PF.plot_recip_min_pair_AS_upset_nt_box(md_data)
#    plt.savefig("test_recip_min_pair_AS_upset_nt_box.png",dpi=600,format="png")
    plt.savefig("{}/recip_min_pair_AS_upset_nt_box.png".format(outdir),dpi=600,format="png")
    plt.clf()
    
    # Upset plot for number of genes with each kind of AS when comparing the
    #   reciprocally min match pairs of the two datasets
    PF.plot_gene_AS_upset(md_data)
#    plt.savefig("test_xcrpt_gene_AS_upset.png",dpi=600,format="png")
    plt.savefig("{}/xcrpt_gene_AS_upset.png".format(outdir),dpi=600,format="png")
    plt.clf()
    
    # Upset plot for number of genes with each kind of AS when comparing the
    #   reciprocally min match pairs of the two datasets with added box plots
    #   of the number of min match pairs present in the gene,
    #   the number of transcripts in each dataset, the avg number of nt different,
    #   and the avg proportion of nt different
    PF.plot_gene_AS_upset_nt_box(md_data,name1,name2)
#    plt.savefig("test_xcrpt_gene_AS_upset_nt_box.png",dpi=600,format="png")
    plt.savefig("{}/xcrpt_gene_AS_upset_nt_box.png".format(outdir),dpi=600,format="png")
    plt.clf()

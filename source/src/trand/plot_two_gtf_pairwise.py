#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:33:04 2021

@author: adalena.nanni
"""

import matplotlib.pyplot as plt

# Import plot functions
from . import plot_functions as PF


def plot_two_gtf_pairwise(outdir, md_data, f1_odds, f2_odds, name1, name2):
    """
    Plot 2 GTF Pairwise output including the following:

        1) Pie charts for the number of genes where the transcripts are the same,
            where one dataset has more transcripts that the onther, and where
            genes are exclusive (pie chart for all genes, single transcript in
            at least one dataset, multi-transcript in at least one dataset, and
            multi-transcript in both datasets)

        2) UpSet plot for number of transcript in genes

        3) A stacked bar chart of the genes with reciprocal minimum pairs for
            all transcript pairs possible (Match:Reciprocal Pairs or Greater:Reciprocal Pairs),
            only some of the transcript pairs (Match:Partial Reciprocal Pairs or Greater:Partial Reciprocal Pairs),
            or none of the transcript pars (Match:No Reciprocal Pairs or Greatre:No Reciprocal Pairs)

        4) The proportional stacked bar chart of (3)

        5) For all minimum pairs of transcripts in either dataset: Upset plot
            for the number of transcripts with each kind of alternative splicing (AS)
            and the distribution of nt differences between the pairs plotted

        6) For minimum pair transcripts of first dataset: Upset plot
            for the number of transcripts with each kind of alternative splicing (AS)
            and the distribution of nt differences between the pairs plotted

        7) For minimum pair transcripts of second dataset: Upset plot
            for the number of transcripts with each kind of alternative splicing (AS)
            and the distribution of nt differences between the pairs plotted
             
        8) For all reciprocal minimum pair transcripts: Upset plot for the number
            of transcripts with each kind of alternative splicing (AS)
            and the distribution of nt differences between the pairs plotted

        9) For all extra minimum pair transcripts: Upset plot for the number
            of transcripts with each kind of alternative splicing (AS)
            and the distribution of nt differences between the pairs plotted
            
        10) Upset plot for number of genes with each kind of AS when comparing the
            reciprocally min match pairs of the two datasets with added box plots
            of the number of min match pairs present in the gene, the number of
            transcripts in each dataset, the average number of nt different,
            and the average proportion of nt different

        11) Plot average number of nt different in recip. min match pairs against
            avg number of nt different in min pairs that are extra
            (either in the greater sets of in set that do not match)

        12) Same as plot (11) with zoom in axises to the mean value of non-zero
            nucleotide differences (max of the two means from recip. min pairs and extras)
    """
# Plot 1
    # Plot pie chart for number of genes where:
    #   the numbers of transcripts are the same (Match) between the two datasets,
    #   where one dataset has more transcripts that the other,
    #   or where one dataset has transcripts and the other has none for the gene

    # Plot pie charts split by 1) all genes, 2) single transcript in at least 1 dataset,
    #   3) multi-transcript in at least 1 dataset, and 4) multi-transcript in both datasets
    PF.plot_transcript_in_gene_split_pie(md_data, f1_odds, f2_odds, name1,
                                         name2, "{}/transcript_in_gene_split_pie.rtf".format(outdir))
    plt.savefig("{}/transcript_in_gene_split_pie.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 2
    # Upset plot for number of transcripts in genes
    PF.plot_transcript_in_gene_upset(md_data, f1_odds, f2_odds, name1, name2,
                                     "{}/transcript_in_gene_upset.rtf".format(outdir))
    plt.savefig("{}/transcript_in_gene_upset.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 3
    # Plot stacked bar chart of genes with reciprocal minimum pairs for:
    #   1) all transcript pairs possible (Match:Reciprocal Pairs or Greater:Reciprocal Pairs)
    #   2) only some of the transcript pairs (Match:Partial Reciprocal Pairs or Greater:Partial Reciprocal Pairs)
    #   3) none of the transcript pars (Match:No Reciprocal Pairs or Greater:No Reciprocal Pairs)
    PF.plot_gene_stack(md_data, name1, name2, "{}/transcript_in_gene_stackCount.rtf".format(outdir))
    plt.savefig("{}/transcript_in_gene_stackCount.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 4
    # With percentage of genes instead of count
    PF.plot_gene_stack(md_data, name1, name2,
                       "{}/transcript_in_gene_stack_proportion.rtf".format(outdir),
                       useProp=True)
    plt.savefig("{}/transcript_in_gene_stack_proportion.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 5
    # For all minimum pairs of transcripts in either dataset:
    # Upset plot for the number of transcripts with each kind of alternative splicing (AS)
    #   and the distribution of nt differences between the pairs plotted
    PF.plot_min_pair_AS_upset_nt_box(md_data, name1, name2,
                                           "{}/all_min_pair_AS_upset_nt_box.rtf".format(outdir),
                                           reciprocal=False, pairs="all")
    plt.savefig("{}/all_min_pair_AS_upset_nt_box.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 6
    # For minimum pair transcripts of first dataset:
    # Upset plot for the number of transcripts with each kind of alternative splicing (AS)
    #   and the distribution of nt differences between the pairs plotted
    PF.plot_min_pair_AS_upset_nt_box(md_data, name1, name2,
                                           "{}/{}_min_pair_AS_upset_nt_box.rtf".format(outdir, name1),
                                           reciprocal=False, pairs="first")
    plt.savefig("{}/{}_min_pair_AS_upset_nt_box.png".format(outdir, name1), dpi=600, format="png")
    plt.clf()

# Plot 7
    # For minimum pair transcripts of second dataset:
    # Upset plot for the number of transcripts with each kind of alternative splicing (AS)
    #   and the distribution of nt differences between the pairs plotted
    PF.plot_min_pair_AS_upset_nt_box(md_data, name1, name2,
                                           "{}/{}_min_pair_AS_upset_nt_box.rtf".format(outdir, name2),
                                           reciprocal=False, pairs="second")
    plt.savefig("{}/{}_min_pair_AS_upset_nt_box.png".format(outdir, name2), dpi=600, format="png")
    plt.clf()

# Plot 8
    # For all reciprocal minimum pair transcripts:
    # Upset plot for the number of transcripts with each kind of AS and the distribution
    #   of nt differences between the pairs plotted above
    PF.plot_min_pair_AS_upset_nt_box(md_data, name1, name2,
                                           "{}/recip_min_pair_AS_upset_nt_box.rtf".format(outdir),
                                           reciprocal=True)
    plt.savefig("{}/recip_min_pair_AS_upset_nt_box.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 9
    # For all extra minimum pair transcripts:
    # Upset plot for the number of transcripts with each kind of AS and the distribution
    #   of nt differences between the pairs plotted above
    PF.plot_min_pair_AS_upset_nt_box(md_data, name1, name2,
                                           "{}/extra_min_pair_AS_upset_nt_box.rtf".format(outdir),
                                           reciprocal=False)
    plt.savefig("{}/extra_min_pair_AS_upset_nt_box.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 10
    # Upset plot for number of genes with each kind of AS when comparing the
    #   reciprocally min match pairs of the two datasets with added box plots
    #   of the number of min match pairs present in the gene,
    #   the number of transcripts in each dataset, the avg number of nt different,
    #   and the avg proportion of nt different
    PF.plot_gene_recip_min_AS_upset(md_data, name1, name2,
                                           "{}/recip_min_gene_AS_upset.rtf".format(outdir))
    plt.savefig("{}/recip_min_gene_AS_upset.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 11
    # Plot average number of nt different in recip. min match pairs against
    #   avg number of nt different in min pairs that are extra
    #   (either in the greater sets of in set that do not match)
    PF.plot_gene_avg_nt_diff_pairs(md_data, name1, name2, "{}/gene_avg_nt_diff_pairs.rtf".format(outdir))
    plt.savefig("{}/gene_avg_nt_diff_pairs.png".format(outdir), dpi=600, format="png")
    plt.clf()

# Plot 12
    # Zoom in axises to the mean value of non-zero nt differences (max of the two means from recip. min pairs and extras)
    PF.plot_gene_avg_nt_diff_pairs(md_data, name1, name2,
                                   "{}/gene_avg_nt_diff_pairs_zoomIn.rtf".format(outdir),
                                   zoomMean=True)
    plt.savefig("{}/gene_avg_nt_diff_pairs_zoomIn.png".format(outdir), dpi=600, format="png")
    plt.clf()

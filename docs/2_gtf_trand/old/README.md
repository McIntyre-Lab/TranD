# _TranD_ 2 GTF Input

![TranD_schematic_B](TranD_schematic_B_white_bg.png)

With two GTF files as input, _TranD_ can be used to calculate and summarize complexity metrics of each GTF file and compare transcript structures across annotations with summarizing visualizations.

<br>

## Complexity
The following complexity metrics are calculated and summarized (minimum, Q1, median, Q3, maximum, mean, and standard deviation) for the single GTF input:
* Transcripts per gene
* Unique exons (exons with unique coordinates) per gene
* Exons per transcripts

<br>

## Pairwise Transcriptome Comparison
Running TranD with two GTF files in __pairwise__ mode (_-e pairwise_, default) produces various plots summarizing the gene and transcript level comparisons of the two annotations.

See examples:
* [_C. elegans_ L1 FLAIR vs. IsoSeq3](celegans_L1_FLAIR_vs_IsoSeq3)
* [_C. elegans_ L1 vs. adult FLAIR](celegans_L1_vs_adult_FLAIR)
* [_Z. mays_ B73 root FLAIR vs. IsoSeq3](maize_B73_root_FLAIR_vs_IsoSeq3)
* [_Z. mays_ B73 root vs. endosperm FLAIR](maize_B73_root_vs_B73_endosperm_FLAIR)
* [_D. melanogaster_ FlyBase r6.17 reference vs. _D. simulans_ FlyBase r2.02 reference (mapped to _D. melanogaster_ genome)](dmel617_vs_dsim202_on_dmel)
* [_D. melanogaster_ FlyBase r6.17 reference vs. _D. simulans_ FlyBase r2.02 reference (mapped to _D. simulans_ genome)](dmel617_vs_dsim202_on_dsim)

<br>

For each gene shared between the two GTF files (require same _gene_id_ value), the structural elements of every unique pair of transcripts across the GTF files are described and quantified to produce various [pairwise distance metrics](../transcript_distance_column_descriptions.xlsx) that are output to a CSV file.
* These metrics allow for nucleotide-level descriptions of genes within each of the combinations of alternative splicing classifications.

<br>

Additionally, pairs are evaluated to determine the minimum pair for each transcript. The [minimum distance variables](../minimum_distance_column_descriptions.xlsx) associated with minimum pairs are written to the same CSV output file. Reciprocal minimum pairs, or pairs of transcripts between the GTF files that are identified as the minimums of each other, are labeled and from these categories, gene-level classifications are also assigned.

<br>

![min_pair_selection](../min_pair_selection.png)

<br>

![recip_min_pair_categories](../recip_min_pair_categories.png)

<br>

![AS_categories](../AS_categories.png)
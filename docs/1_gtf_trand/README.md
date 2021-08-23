
# _TranD_ 1 GTF Input

![TranD_schematic_A](TranD_schematic_A_white_bg.png)

With a single GTF as input, _TranD_ can be used to generate a splice-match consolidated transcriptome annotation, calculate and summarize complexity metrics, summarize structural elements of the transcriptome annotation, and summarize transcriptome-wide alternative splicing.

<br>

## Splice-Match Consolidation

To focus on the comparison of structural variation in transcript splicing, transcriptome anntoations can be splice-match consolidated (i.e., transcripts with identical splice junctions are consolidated to one uniquely spliced representative transcript that spans the exon space of the contributing transcripts).

<br>

## Complexity
The following complexity metrics are calculated and summarized (minimum, Q1, median, Q3, maximum, mean, and standard deviation) for the single GTF input:
* Transcripts per gene
* Unique exons (exons with unique coordinates) per gene
* Exons per transcripts

<br>

## Gene Mode
Running TranD with a single GTF input in __gene__ mode (_-e gene_) produces a figure for summarizing the structural elements of the annotation.

See examples:

* [_C. elegans_ WBcel235 Reference Transcriptome Annotation TranD Gene](celegans_WBcel235_ref_trand_gene)
* [_Z. mays_ Mo17 YAN Reference Transcriptome Annotation TranD Gene](maize_mo17_YAN_ref_trand_gene)
* [_C. elegans_ WBcel235 Splice-Match Consolidated Reference Transcriptome Annotation TranD Gene](celegans_WBcel235_splice_match_consol_ref_trand_gene)
* [_Z. mays_ Mo17 YAN Splice-Match Consolidated Reference Transcriptome Annotation TranD Gene](maize_mo17_YAN_splice_match_consol_ref_trand_gene)

<br>

## Pairwise Mode
Running TranD with one GTF input in __pairwise__ mode (_-e pairwise_, default) produces a figure for summarizing the combinations of different alternative splicing events (see below) within multi-transcript genes of the given annotation.

<br>

![AS_categories](../AS_categories.png)

<br>

See examples:

* [_C. elegans_ WBcel235 Reference Transcriptome Annotation TranD Pairwise](celegans_WBcel235_ref_trand_pairwise)
* [_Z. mays_ Mo17 YAN Reference Transcriptome Annotation TranD Pairwise](maize_mo17_YAN_ref_trand_pairwise)
* [_C. elegans_ WBcel235 Splice-Match Consolidated Reference Transcriptome Annotation TranD Pairwise](celegans_WBcel235_splice_match_consol_ref_trand_pairwise)
* [_Z. mays_ Mo17 YAN Splice-Match Consolidated Reference Transcriptome Annotation TranD Pairwise](maize_mo17_YAN_splice_match_consol_ref_trand_pairwise)

<br>

For each gene, the structural elements of each unique pair of transcripts are described and quantified to produce various [pairwise distance metrics](../transcript_distance_column_descriptions.xlsx) that are output to a CSV file.

* These metrics allow for nucleotide-level descriptions of genes within each of the combinations of alternative splicing classifications.

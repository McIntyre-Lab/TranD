# TranD: Transcript Distances

TranD is a collection of tools to facilitate metrics of structural variation for whole genome
transcript annotation files (GTF) that pinpoint structural variation to the nucleotide level.

TranD (Transcript Distances) can be used to calculate metrics of structural variation within and
between annotation files (GTF). Structural variation reflects organismal complexity and three
summary metrics for genome level complexity are calculated for every gene in a GTF file:  1) the
number of transcripts per gene; 2) the number of exons per transcript; and 3) the number of
unique exons (exons with unique genomic coordinates) per gene. From each these metrics
distributions a summary statistics such as mean, median, variance, inter-quartile range are
calculated. With 1GTF file input, gene mode can be used to generates these metrics for each gene
and summary statistics and distributions across genes. Distributions are visualized in a series
of plots.  For 1 GTF and 2GTF a pairwise mode calculates distance metrics between 2 transcripts
to the nucleotide.  In 1 GTF this is all possible pairs within the gene and in 2 GTF model this
is all possible pairs among GTF files. The distribution of these metrics across genes are
visualized and summary statistics for structural variations between pairs are calculated and
reported.  Visualizations of the distributions of the frequency of intron retention, alternative
exon usage, donor/acceptor variation and 5', 3' variation in UTR regions are provided as well as
tabular formatted nucleotide level distances for both 1GTF and 2 GTF.

## Installation

To install TranD run
```
pip install trand
```
in activated conda or virtualenv environment

If you clone this repository you can use the source/requirements.txt file to pip install required dependencies with
```
pip install -r source/requirements.txt
```
See the [PyPI TranD page](https://pypi.org/project/trand/) for reference.

We will add a conda trand package in the near future.

# TranD: Transcript Distances

**Below is a basic description of the functionality of TranD. New users, please check the [Wiki](https://github.com/McIntyre-Lab/TranD/wiki). The [User Guide](https://github.com/McIntyre-Lab/TranD/wiki/User-Guide) contains more detailed information on how to run and install TranD.**

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

The easiest way to install TranD is by using conda or mamba.

```
conda install -c bioconda trand
```

If you already have a suitable environment with bedtools installed and only want to add python
packages run
```
pip install trand
```
in an activated conda or virtualenv environment

If you clone this repository you can use the source/requirements.txt file to pip install required dependencies with
```
pip install -r source/requirements.txt
```
See the [PyPI TranD page](https://pypi.org/project/trand/) for reference.

We will add a conda trand package in the near future.

## Folder Description

- **TMM_input** and **TMM_output**: Input and Output used in the [Transcript Model Map](https://github.com/McIntyre-Lab/TranD/wiki/Create-within-Species-Transcript-Model-Map) example for the wiki.
- **docs**: Documentation for with details on the process of creating the [Pre-Computed Files](https://github.com/McIntyre-Lab/TranD/wiki/Precomputed-Files)
- **source**: Contains the python code for TranD
- **tests**: Contains shell scripts used to test TranD during development
- **utilities**: Contains the python code for the various utilities to be used in conjunction with TranD. For more information please see the [Utility Descriptions](https://github.com/McIntyre-Lab/TranD/wiki/Utility-Descriptions-(with-Examples))
- **utility_input** and **utility_output**: Input and Output used as examples in the [Utility Descriptions](https://github.com/McIntyre-Lab/TranD/wiki/Utility-Descriptions-(with-Examples))

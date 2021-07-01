#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA=data/test/test_exons_gene.gtf

python src/eventanalysis/event_analysis.py -f -e gene -o testout/gene_small_synthetic "$@" ${DATA}

date

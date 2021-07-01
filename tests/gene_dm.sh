#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA=data/raw/mel_head_PB.pb2fbgn.sorted.gtf

time python src/eventanalysis/event_analysis.py -e gene -o testout/gene_full_mel_head "$@" ${DATA}

date

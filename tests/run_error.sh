#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

#Test1: Too many arguments
python eventanalysis/src/event_analysis.py \
    data/raw/mel_head_PB.pb2fbgn.sorted.gtf \
    data/raw/mel_head_PB.pb2fbgn.sorted.duplicate.copy.gtf \
    data/raw/mel_head_PB.pb2fbgn.sorted.triplicate.copy.gtf

# TODO
#Test1: fasta data in a file (when using gff)

date

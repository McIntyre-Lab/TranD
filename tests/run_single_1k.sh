#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

python eventanalysis/src/event_analysis.py -o testout "$@" data/test/mel_head_PB.pb2fbgn.sorted.1000.gtf

date

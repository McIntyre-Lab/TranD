#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA=data/test/test_two_identical.gtf

python eventanalysis/src/event_analysis.py -o testout "$@" ${DATA}

date

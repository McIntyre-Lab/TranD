#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA="data/raw/test_set_B73_END.gtf data/raw/test_set_B73_R.gtf"

python src/eventanalysis/event_analysis.py -f -o testout "$@" ${DATA}

date

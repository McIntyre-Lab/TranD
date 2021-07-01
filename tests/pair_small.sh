#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA=data/test/test_exons.gtf

CMD="python src/eventanalysis/event_analysis.py -o testout/pair_small $@ ${DATA}"
echo "CMD: ${CMD}"
eval ${CMD}

date

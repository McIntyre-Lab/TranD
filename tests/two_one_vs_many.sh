#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/2file_one_vs_many_exons_1.gtf data/test/2file_one_vs_many_exons_2.gtf"

CMD="python src/eventanalysis/event_analysis.py -f -o testout/two_files_one_vs_many $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/identical_1_FLAIR_R_Zm00001d037333_2_d1.gtf data/test/identical_2_FLAIR_END_Zm00001d037333_3_d2.gtf"

CMD="python src/eventanalysis/event_analysis.py -f -o testout/two_files_identical $@ ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

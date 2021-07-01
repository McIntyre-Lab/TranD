#!/bin/bash

export PATH=${HOME}/projects/bioinformatics/mcintyre-cid/Event_Analysis_2.0/conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/2file_missing_fragments_B73_R_FLAIR.gtf data/test/2file_missing_fragments_B73_END_FLAIR.gtf"

CMD="python src/eventanalysis/event_analysis.py -f -o testout/two_files_missing_fragments $@ ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

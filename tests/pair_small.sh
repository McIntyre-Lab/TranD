#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/test_exons.gtf"

CMD="python src/eventanalysis/event_analysis.py -o testout/pair_small $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

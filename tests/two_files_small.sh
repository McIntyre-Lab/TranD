#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/raw/test_set_B73_END.gtf data/raw/test_set_B73_R.gtf"

CMD="trand -f -o testout/two_files_small $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date
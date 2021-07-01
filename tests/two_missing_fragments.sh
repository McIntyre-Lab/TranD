#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/2file_missing_fragments_B73_R_FLAIR.gtf data/test/2file_missing_fragments_B73_END_FLAIR.gtf"

CMD="trand -f -o testout/two_files_missing_fragments $@ ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

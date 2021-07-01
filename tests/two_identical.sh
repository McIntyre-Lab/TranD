#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/identical_1_FLAIR_R_Zm00001d037333_2_d1.gtf data/test/identical_2_FLAIR_END_Zm00001d037333_3_d2.gtf"

CMD="trand -f -o testout/two_files_identical $@ ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

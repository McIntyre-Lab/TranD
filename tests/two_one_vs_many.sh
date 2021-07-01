#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/2file_one_vs_many_exons_1.gtf data/test/2file_one_vs_many_exons_2.gtf"

CMD="trand -f -o testout/two_files_one_vs_many $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

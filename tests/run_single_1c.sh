#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/mel_head_PB.pb2fbgn.sorted.100.gtf"

CMD="trand -o testout/single_1c $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

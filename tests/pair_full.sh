#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA=data/raw/mel_head_PB.pb2fbgn.sorted.gtf

TEST="pair_full"

CMD="trand -f -o testout/${TEST} $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

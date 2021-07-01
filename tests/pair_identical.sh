#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/test_two_identical.gtf"
TEST="pair_identical"

CMD="trand -f -o testout/${TEST} $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

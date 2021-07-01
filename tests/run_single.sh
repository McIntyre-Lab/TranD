#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

CMD="trand -f -o testout/run_single $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

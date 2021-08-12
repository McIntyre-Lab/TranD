#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

TEST="pair_small"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

DATA="test_exons.gtf"

CMD="trand -o testout/pair_small $* ${INPUT_DIR}/${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date
echo "### Finished test: ${TEST} on $(date)"

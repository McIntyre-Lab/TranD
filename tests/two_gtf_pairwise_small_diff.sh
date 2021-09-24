#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

TEST="two_gtf_pairwise_small"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
echo "### Starting test: ${TEST}"

DATA="${INPUT_DIR}/test_set_B73_END.gtf ${INPUT_DIR}/test_set_B73_R_diff.gtf"

CMD="trand -o ${OUTPUT_DIR} $* ${DATA}"
echo "CMD: ${CMD}"
time eval "${CMD}"

date
echo "### Finished test: ${TEST} on $(date)"

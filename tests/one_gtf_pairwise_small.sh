#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

TEST="one_gtf_pairwise_small"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

DATA="test_exons.gtf"

#Performance profiling
EXE="${HOME}/projects/bioinformatics/mcintyre-cid/TranD/conda/bin/trand"
CMD="python -m cProfile -o ${OUTPUT_DIR}/cprofile.log ${EXE} -f -s -o ${OUTPUT_DIR} $* ${INPUT_DIR}/${DATA}"
#standard test
#CMD="trand -f -s -o ${OUTPUT_DIR} $* ${INPUT_DIR}/${DATA}"
echo "CMD: ${CMD}"
time eval "${CMD}"

date
echo "### Finished test: ${TEST} on $(date)"

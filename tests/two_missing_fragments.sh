#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

TEST="gene_small_pairwise"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

DATA="${INPUT_DIR}/2file_missing_fragments_B73_R_FLAIR.gtf ${INPUT_DIR}/2file_missing_fragments_B73_END_FLAIR.gtf"

CMD="trand -f -o testout/two_files_missing_fragments $@ ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

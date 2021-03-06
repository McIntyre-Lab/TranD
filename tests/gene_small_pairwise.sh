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

DATA="test_exons_gene.gtf"

CMD="trand -f -e pairwise -s -o ${OUTPUT_DIR} $* ${INPUT_DIR}/${DATA}"
echo "CMD: ${CMD}"
time eval "${CMD}"

date

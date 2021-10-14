#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

TEST="gene_small_full_gene"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

DATA1="identical_1_FLAIR_R_Zm00001d037333_2_d1.gtf"
DATA2="identical_2_FLAIR_END_Zm00001d037333_3_d2.gtf"

CMD="trand -f -s -o testout/two_files_identical $@ ${INPUT_DIR}/${DATA1} ${INPUT_DIR}/${DATA2}"
echo "CMD: ${CMD}"
eval time "${CMD}"

date

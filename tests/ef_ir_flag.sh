#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

TEST="one_gtf_ef_ir_flags_freq"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

#DATA="frag_IR_and_xcrpt_assignment_bug_ex.gtf"
DATA="Zm00009a000090.gtf"

CMD="trand -o ${OUTPUT_DIR} -f -k -e gene -s $* ${INPUT_DIR}/${DATA}"
echo "CMD: ${CMD}"
time eval "${CMD}"

date
echo "### Finished test: ${TEST} on $(date)"


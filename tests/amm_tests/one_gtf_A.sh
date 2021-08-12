#!/bin/bash

## A) calculate complexity measures only (do not generate events) on a GTF file without consolidation, overwrite any existing output files:

PROJ=/TB14/TB14/github/TranDi_EA
SCRIPT=$PROJ/source/src

export PATH=/TB14/TB14/conda_envs/tranD/bin:${PATH}

pwd;hostname;date

INPUT=$PROJ/galaxy/test-data
ONE_GTF=cel_adult_FLAIR_chr4

OUTPUT=$PROJ/galaxy/amm_testout
    mkdir -p $OUTPUT

TEST="one_gtf_A"
TESTOUT=$OUTPUT/${TEST}

rm -rf ${TESTOUT}
mkdir -p ${TESTOUT}

echo "### Starting test: ${TEST}"

trand $INPUT/${ONE_GTF}.gtf --complexityOnly --outdir $TESTOUT --force

date
echo "### Finished test: ${TEST} on $(date)"

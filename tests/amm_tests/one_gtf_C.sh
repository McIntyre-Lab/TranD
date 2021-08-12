#!/bin/bash


## C) calculate complexity measures on consolidated GTF file, generate event file by comparing all transcripts for a gene and create a transcriptome summary plot 
    ## (note that if you want to keep intron retention events add the --keepir):

PROJ=/TB14/TB14/github/TranDi_EA
SCRIPT=$PROJ/source/src

export PATH=/TB14/TB14/conda_envs/tranD/bin:${PATH}

pwd;hostname;date

INPUT=$PROJ/galaxy/test-data
ONE_GTF=cel_adult_FLAIR_chr4

OUTPUT=$PROJ/galaxy/amm_testout
    mkdir -p $OUTPUT

TEST="one_gtf_C"
TESTOUT=$OUTPUT/${TEST}

rm -rf ${TESTOUT}
mkdir -p ${TESTOUT}

echo "### Starting test: ${TEST}"

trand $INPUT/${ONE_GTF}.gtf --consolidate --consolPrefix tr --ea gene --outdir $TESTOUT --force --logfile $TESTOUT/${TEST}.log

date
echo "### Finished test: ${TEST} on $(date)"

#!/bin/bash



## D) calculate complexity measures on consolidated GTF file, generate event file by a pairwise comparison of transcripts for a gene
    ## (keeping intron retention events) and calculate distances

PROJ=/TB14/TB14/github/TranDi_EA
SCRIPT=$PROJ/source/src

export PATH=/TB14/TB14/conda_envs/tranD/bin:${PATH}

pwd;hostname;date

INPUT=$PROJ/galaxy/test-data
ONE_GTF=cel_adult_FLAIR_chr4

OUTPUT=$PROJ/galaxy/amm_testout
    mkdir -p $OUTPUT

TEST="one_gtf_D"
TESTOUT=$OUTPUT/${TEST}

rm -rf ${TESTOUT}
mkdir -p ${TESTOUT}

echo "### Starting test: ${TEST}"

trand $INPUT/${ONE_GTF}.gtf --consolidate --consolPrefix tr --ea pairwise --keepir --outdir $TESTOUT --force --logfile $TESTOUT/${TEST}.log

date
echo "### Finished test: ${TEST} on $(date)"

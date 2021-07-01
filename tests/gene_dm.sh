#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA=data/raw/mel_head_PB.pb2fbgn.sorted.gtf

trand -e gene -o testout/gene_full_mel_head "$@" ${DATA}

date

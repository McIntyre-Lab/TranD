#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

#Test1: Too many arguments
trand \
    data/raw/mel_head_PB.pb2fbgn.sorted.gtf \
    data/raw/mel_head_PB.pb2fbgn.sorted.duplicate.copy.gtf \
    data/raw/mel_head_PB.pb2fbgn.sorted.triplicate.copy.gtf

# TODO
#Test1: fasta data in a file (when using gff)

date

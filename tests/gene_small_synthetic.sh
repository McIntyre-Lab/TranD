#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=${SCRIPT_DIR}/../../conda/bin:${PATH}
pwd;hostname;date

DATA="data/test/test_exons_gene.gtf"

TEST="gene_small_synthetic"

CMD="trand -f -e gene -o testout/${TEST} $* ${DATA}"
echo "CMD: ${CMD}"
eval "${CMD}"

date

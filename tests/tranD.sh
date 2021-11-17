#!/bin/bash
#    <test>
#       <!-- Calculate complexity measures only (do not generate events) on a GTF file after consolidation -->
#       <param name="gtf1Dataset" value="cel_adult_FLAIR_chr4.gtf"/>
#       <param name="complexityOnly" value="True" />
#       <param name="keepir" value="True" />
#       <output_collection name="counts" type="list">
#         <element name="transcriptome_complexity_counts" file="transcriptome_complexity_counts.csv" ftype="csv" />
#       </output_collection> 
#       <output_collection name="plots" type="list">
#         <element name="complexity_plots" file="complexity_plots.png" ftype="png"/>
#       </output_collection>
#       <output_collection name="figureLegend" type="list">
#         <element name="complexity_fig_legend" file="complexity_plots.rtf" ftype="rtf"/>
#       </output_collection>
#    </test>

SCRIPT=$(basename "${BASH_SOURCE[0]}");
TEST="${SCRIPT%.*}"
echo "Echoing ${BASH_SOURCE[0]}"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"
if [[ $# -gt 0 ]]; then OUTPUT_DIR=$1 ; fi
mkdir -p "${OUTPUT_DIR}"

main.py \
    $INPUT_DIR/cel_adult_FLAIR_chr4.gtf \
    --complexityOnly \
    --keepir \
    --outdir $OUTPUT_DIR

echo "### Finished test: ${TEST} on $(date)"

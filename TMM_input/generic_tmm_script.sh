#!/bin/bash

# Below is a code version of the wiki page describing each step of the creation of the
# Transcript Model Map. The variables at the very beginning are provided for
# ease when analyzing any two GTF files.

# For any confusion on any step below, please check the wiki
# for detailed explanation.


SCRIPTS=/location/of/TranD/utilities

GTF1=/path/to/GTF1
GTF2=/path/to/GTF2

NAME1=gtf1
NAME2=gtf2

OUTD=/path/to/output

NUM_CPUS=8

# Calculation for this threshold can be found in the wiki
NT_THRESHOLD=15


echo "Identifying UJCs..."
# ID UJCs
python ${SCRIPTS}/id_ujc.py \
-g ${GTF1} \
-o ${OUTD} \
-p ${NAME1}


python ${SCRIPTS}/id_ujc.py \
-g ${GTF1} \
-o ${OUTD} \
-p ${NAME2}


echo "Creating Union Reference..."
# Create Union Reference
cat ${OUTD}/${NAME1}_UJC.gtf \
    ${OUTD}/${NAME2}_UJC.gtf \
    > ${OUTD}/${NAME1}_${NAME2}_UJC.gtf

python ${SCRIPTS}/id_ujc.py \
-g ${NAME1}_${NAME2}_UJC.gtf \
-o ${OUTD} \
-p ${NAME1}_${NAME2}_union


echo "Running TranD 2GTF for initial pairwise distance calculations..."
# Minimum Pairwise Distance Calculations
trand \
${OUTD}/${NAME1}_UJC.gtf \
${OUTD}/${NAME2}_UJC.gtf \
-o ${OUTD}/TranD_2GTF_output \
-1 ${NAME1} \
-2 ${NAME2} \
-e pairwise \
-p both \
-f \
-n ${NUM_CPUS}

# New Folder: TranD Output
TD_OUT=${OUTD}/TranD_2GTF_output

echo "Running pair classification..."
# Pair Classification
python ${SCRIPTS}/pair_classfication.py \
-d ${TD_OUT}/${NAME1}_vs_${NAME2}_minimum_pairwise_distance.csv \
-s ${NT_THRESHOLD}
-o ${OUTD}/${NAME1}_${NAME2}_pair_classification.csv

# Pairwise Distance for Genes with 1 GTF (only done if the files are not empty)
# Condition runs trand if file is not empty

echo "Checking GTF1 only and GTF2 only..."
if [ -s "${OUTD}/TranD_2GTF_output/${NAME1}_vs_${NAME2}_gtf1_only.gtf" ]; then

	trand \
	${TD_OUT}/${NAME1}_vs_${NAME2}_gtf1_only.gtf \
	-o ${OUTD}/TranD_GTF1_only_output \
	-e pairwise \
	-f \
	-n ${NUM_CPUS}

else
	echo "The GTF1 only file is empty."

fi

if [ -s "${OUTD}/TranD_2GTF_output/${NAME1}_vs_${NAME2}_gtf2_only.gtf" ]; then

	trand \
	${TD_OUT}/${NAME1}_vs_${NAME2}_gtf2_only.gtf \
	-o ${OUTD}/TranD_GTF2_only_output \
	-e pairwise \
	-f \
	-n ${NUM_CPUS}

else
	echo "The GTF2 only file is empty."

fi

echo "Running ID_ERG..."
# Run ID_ERG

python ${SCRIPTS}/id_ERG.py \
-i ${TD_OUT}/${NAME1}_vs_${NAME2}_minimum_pairwise_distance.csv \
-g \
-p ${NAME1}_vs_${NAME2} \
-o ${OUTD}/ERG/2GTF \
-ir Y

# Run ID_ERG on GTF1 only and GTF2 only, if necessary
if [ -s "${OUTD}/TranD_2GTF_output/${NAME1}_vs_${NAME2}_gtf1_only.gtf" ]; then

	python ${SCRIPTS}/id_ERG.py \
	-i ${OUTD}/TranD_GTF1_only_output/pairwise_distance.csv \
	-g \
	-p ${NAME1}_vs_${NAME2}_gtf1_only \
	-o ${OUTD}/ERG/GTF1_only \
	-ir Y \
	-w 1

else
	echo "The GTF1 only file is empty. No ERGs identified"

fi

if [ -s "${OUTD}/TranD_2GTF_output/${NAME1}_vs_${NAME2}_gtf2_only.gtf" ]; then

	python ${SCRIPTS}/id_ERG.py \
	-i ${OUTD}/TranD_GTF2_only_output/pairwise_distance.csv \
	-g \
	-p ${NAME1}_vs_${NAME2}_gtf2_only \
	-o ${OUTD}/ERG/GTF2_only \
	-ir Y \
	-w 2
else
	echo "The GTF2 only file is empty. No ERGs identified"

fi

# ERG Filepaths:
ERG_2GTF=${OUTD}/ERG/2GTF
ERG_GTF1=${OUTD}/ERG/GTF1_only
ERG_GTF2=${OUTD}/ERG/GTF2_only

echo "Concatenating ERG fie
# Create Union ERG Files

cat ${ERG_2GTF}/${NAME1}_vs_${NAME2}_erg_gtf.gtf \
    ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf1_only_erg_gtf.gtf \
    ${ERG_GTF2}/${NAME1}_vs_${NAME2}_gtf2_only_erg_gtf.gtf \
    2>/dev/null \ # Allows for concatenation even if one file does not exist
    > ${OUTD}/ERG/${NAME1}_vs_${NAME2}_union_erg_gtf.gtf


cat ${ERG_2GTF}/${NAME1}_vs_${NAME2}_xscript_output.csv \
    <(tail +2 ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf1_only_xscript_output.csv) \
    <(tail +2 ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf2_only_xscript_output.csv) \
    2>/dev/null \ # Allows for concatenation even if one file does not exist
    > ${OUTD}/ERG/${NAME1}_vs_${NAME2}_union_xscript_output.csv


cat ${ERG_2GTF}/${NAME1}_vs_${NAME2}_gene_output.csv \
    <(tail +2 ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf1_only_gene_output.csv) \
    <(tail +2 ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf2_only_gene_output.csv) \
    2>/dev/null \ # Allows for concatenation even if one file does not exist
    > ${OUTD}/ERG/${NAME1}_vs_${NAME2}_union_gene_output.csv


cat ${ERG_2GTF}/${NAME1}_vs_${NAME2}_erg_output.csv \
    <(tail +2 ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf1_only_erg_output.csv) \
    <(tail +2 ${ERG_GTF1}/${NAME1}_vs_${NAME2}_gtf2_only_erg_output.csv) \
    2>/dev/null \ # Allows for concatenation even if one file does not exist
    > ${OUTD}/ERG/${NAME1}_vs_${NAME2}_union_erg_output.csv

echo "Transcript Model Map Complete!"

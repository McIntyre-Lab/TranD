#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:29:04 2024

@author: k.bankole
"""

import argparse
import pandas as pd
import trand.io
import time
import numpy as np
import os


def getOptions():
    """

    Function to store user input via argparse

    Returns
    -------
    args : ARGPARSE ARGUMENTS
            User input via argparse.

    """

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="A script that compares a desired GTF (-d) (typically "
                                     "reads in GTF form) to an exon segment (ES) "
                                     "GTF (-er). Creates exon segment patterns (ESP) which are "
                                     "binary patterns indicating which of a gene's exon "
                                     "segments a transcript has exons within. Outputs two files: "
                                     "to desired output directory (-o). "
                                     "1. A list of transcripts and their ESPs. 2. A flag file "
                                     "that indicates which of the gene's exon segments the "
                                     "transcript has. Script will also output a list of genes"
                                     "that only appear in one GTF. There is also an option to add an "
                                     "output prefix (-x).")

    # INPUT
    parser.add_argument(
        "-es",
        "--es-gtf",
        dest="esFile",
        required=True,
        help="Location of ES GTF"
    )

    parser.add_argument(
        "-d",
        "--data-gtf",
        dest="dataFile",
        required=True,
        help="Location of data GTF"
    )

    # OUTPUT
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output file path. Must already exist."
    )

    parser.add_argument(
        "-x",
        "--prefix",
        dest="prefix",
        required=False,
        help="Output prefix."
    )

    parser.add_argument(
        "-s",
        "--sample-ID",
        dest="sampleID",
        required=False,
        help="Optional SampleID. Will create a sampleID column in output."
    )

    args = parser.parse_args()
    return args


def main():

    # NOTE: This script has the **exact same logic** as the ERP version.
    esFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc_es.gtf"
    dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_noMultiGene.gtf"

    outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # prefix = "test_prefix"
    sampleID = None
    prefix = None

    esFile = args.esFile
    dataFile = args.dataFile
    outdir = args.outdir
    prefix = args.prefix
    sampleID = args.sampleID

    alphatic = time.perf_counter()

    # Read in both GTFs and subset them to genes that are in both GTF
    inGeneDf = trand.io.read_exon_data_from_file(esFile)
    inDataDf = trand.io.read_exon_data_from_file(dataFile)

    uniqDataGeneSet = set(inDataDf['gene_id'])
    uniqRefGeneSet = set(inGeneDf['gene_id'])

    # Store genes only in one GTF for later output
    refOnlyGnLst = list(uniqRefGeneSet - uniqDataGeneSet)
    dataOnlyGnLst = list(uniqDataGeneSet - uniqRefGeneSet)

    genesInBoth = list(uniqRefGeneSet.intersection(uniqDataGeneSet))

    geneDf = inGeneDf[inGeneDf['gene_id'].isin(genesInBoth)].copy()
    dataDf = inDataDf[inDataDf['gene_id'].isin(genesInBoth)].copy()

    # geneDf = inGeneDf[inGeneDf['gene_id'].isin(["LOC120456871"])].copy()
    # dataDf = inDataDf[inDataDf['gene_id'].isin(["LOC120456871"])].copy()

    # Clean up ES GTF dataframe (geneDf)
    geneDf = geneDf[['gene_id', 'seqname', 'start', 'end', 'strand']].copy()
    geneDf = geneDf.sort_values(
        ['seqname', 'gene_id', 'start'], ignore_index=True)

    # Check that each gene is only on one strand (don't know why they wouldn't be)
    singleStrandGene = geneDf.groupby('gene_id').agg(
        set)['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandGene.all():
        print("There are genes belonging to more than one strand. Quitting.")
        quit()

    # Assign each exon in the ES GTF its ES ID
    geneDf['ES'] = geneDf['gene_id'] + ':ES' + \
        (geneDf.groupby('gene_id').cumcount() + 1).astype(str)

    # Create a dictionary of genes and their ESs. Sort ESIDs to be in numerical order (matches 5'->3' relative to + strand)
    geneDct = dict(geneDf.groupby('gene_id').apply(
        lambda x: sorted(set(x['ES']), key=lambda x: int(x.split("ES")[1]))))

    # TODO: CHECK THAT ALL SETS ARE OF SIZE ONE
    # Create a dictionary of ESs and their information
    esDf = geneDf.groupby('ES').agg('first')
    esDf['length'] = esDf['end'] - esDf['start']
    esDct = esDf.to_dict(orient='index')

    dataDf['numExon'] = dataDf.groupby('transcript_id')[
        'transcript_id'].transform('count')

    dataDf['dataOnlyExon'] = np.nan

    records = dataDf.to_dict('records')

    for row in records:

        gene = row['gene_id']
        # jxnHash = row['transcript_id']

        matchingESIDLst = []

        # if gene in geneDct.keys():
        for esID in geneDct.get(gene):
            # print(esID)
            esInfo = esDct.get(esID)
            # print(esInfo)

            # print("looping...")

            if max(row['start'], esInfo['start']) < min(row['end'], esInfo['end']):
                # print(row)
                # print(esID)
                # print(erInfo)

                matchingESIDLst.append(esID)

        if matchingESIDLst:
            row['ES'] = matchingESIDLst
        else:
            row['dataOnlyExon'] = "{}:{}_{}".format(
                gene, row['start'], row['end'])

    dataWithESDf = pd.DataFrame(records)

    intmdDf = dataWithESDf[['seqname', 'gene_id', 'transcript_id', 'ES',
                            'dataOnlyExon', 'numExon', 'strand']]
    intmdDf = intmdDf.explode('ES')

    xscriptESDf = intmdDf.groupby('transcript_id').agg({
        'ES': lambda x: set(x.dropna()),
        'numExon': max,
        'dataOnlyExon': lambda x: set(x.dropna()),
        'strand': set,
        'seqname': set
    }).reset_index()

    # Check for no multi-strand transcripts
    singleStrandXscript = xscriptESDf['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandXscript.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        xscriptESDf['strand'] = xscriptESDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrXscript = xscriptESDf['seqname'].apply(lambda x: len(x) == 1)

    if not singleChrXscript.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        xscriptESDf['seqname'] = xscriptESDf['seqname'].apply(
            lambda x: list(x)[0])

    xscriptESDct = dict(zip(xscriptESDf['transcript_id'], xscriptESDf['ES']))

    dataOnlyExonDct = dict(
        zip(xscriptESDf['transcript_id'], xscriptESDf['dataOnlyExon']))

    # dataWithESDf = dataWithESDf[
    #     (dataWithESDf['transcript_id'] ==
    #       "34a6dd0389208ff91783b8ec557da2427785a0988ffa8243635e75c13b7f334c")
    #     | (dataWithESDf['transcript_id'] == "068509f23790d261052860383c257e94c8d917c9c67a0140875242cd7302b738")]

    loopLst = [tuple(x) for x in dataWithESDf[[
        'gene_id', 'transcript_id', 'strand', 'seqname']].drop_duplicates().to_records(index=False)]

    xscriptLst = []
    geneLst = []
    seqnameLst = []
    erLst = []
    flagLst = []
    strandLst = []
    lngthLst = []

    patternDct = dict()
    for gene, transcript, strand, seqname in loopLst:

        geneESLst = geneDct.get(gene)
        xscriptESSet = xscriptESDct.get(transcript)

        pttrnLst = ["1" if ES in xscriptESSet else "0" for ES in geneESLst]
        esIDLst = [ES for ES in geneESLst if ES in xscriptESSet]

        if strand == "-":
            pttrnLst.reverse()
            esIDLst.reverse()

        pattern = strand + "_" + ''.join(map(str, pttrnLst))
        patternESID = strand + "_" + "_".join(esIDLst)

        patternDct[transcript] = [pattern, patternESID, gene]

        for exonRegion in geneESLst:

            if exonRegion in xscriptESSet:
                flag = 1
            else:
                flag = 0

            xscriptLst.append(transcript)
            geneLst.append(gene)
            seqnameLst.append(seqname)
            erLst.append(exonRegion)
            flagLst.append(flag)
            strandLst.append(strand)
            lngthLst.append(esDct[exonRegion]['length'])

        if dataOnlyExonDct[transcript]:

            for exonRegion in dataOnlyExonDct[transcript]:
                xscriptLst.append(transcript)
                geneLst.append(gene)
                seqnameLst.append(seqname)
                erLst.append(exonRegion)
                flagLst.append(1)
                strandLst.append(strand)

                # TODO: CHANGE THIS IF THE DATA ONLY EXON FORMAT CHANGES
                startNEnd = exonRegion.split(':')[1].split('_')
                exonLngth = int(startNEnd[1]) - int(startNEnd[0])

                lngthLst.append(exonLngth)

    outFlagDf = pd.DataFrame({
        'jxnHash': xscriptLst,
        'geneID': geneLst,
        'seqname': seqnameLst,
        'strand': strandLst,
        'ES': erLst,
        'flagES': flagLst,
        'lengthES': lngthLst
    })

    # Making pattern output file
    pttrnInfo = [(xscript, *info) for xscript, info in patternDct.items()]
    patternDf = pd.DataFrame(pttrnInfo, columns=[
        'transcript_id', 'ESP', 'patternESID', 'geneID'])

    outPatternDf = pd.merge(xscriptESDf, patternDf, on=[
        'transcript_id'], how='outer', indicator='merge_check')
    outPatternDf.rename(columns={'transcript_id': 'jxnHash'}, inplace=True)

    if not (outPatternDf['merge_check'] == 'both').all():
        raise Exception(
            "Something went wrong. Merge of patterns and xscript information failed.")

    outPatternDf['flagDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: len(x) != 0)

    # TODO: no numES in ESP file
    outPatternDf['numDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(len)

    outPatternDf['dataOnlyESID'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    outFlagDf = outFlagDf.sort_values(by=['geneID', 'jxnHash'])
    outPatternDf = outPatternDf.sort_values(by=['geneID', 'jxnHash'])

    outPatternDf = outPatternDf[['jxnHash', 'geneID', 'seqname', 'strand',
                                 'ESP', 'patternESID', 'numExon',
                                 'numDataOnlyExon', 'dataOnlyESID']]

    if sampleID:
        outFlagDf['sampleID'] = sampleID
        outPatternDf['sampleID'] = sampleID

    # Do not uncomment. Will probably crash the script.
    # wideDf = pd.pivot_table(outDf, values='flagES', index=['jxnHash','geneID'], columns='exonRegion', fill_value=0)

    # Output
    esName = os.path.splitext(os.path.basename(esFile))[0]
    dataName = os.path.splitext(os.path.basename(dataFile))[0]

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    erpFile = outPrefix + "{}_vs_{}_ESP.csv".format(esName, dataName)

    outPatternDf.to_csv(
        erpFile, index=False)

    flagFile = outPrefix + "{}_vs_{}_flagES.csv".format(esName, dataName)

    outFlagDf.to_csv(
        flagFile, index=False)

    if refOnlyGnLst:
        pd.Series(refOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_anno_only_genes.txt".format(esName, dataName), index=False, header=False)

    if dataOnlyGnLst:
        pd.Series(dataOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_data_only_genes.txt".format(esName, dataName), index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

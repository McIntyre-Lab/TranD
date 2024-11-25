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
import re


def getOptions():
    """

    Function to store user input via argparse

    Returns
    -------
    args : ARGPARSE ARGUMENTS
            User input via argparse.

    """

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="A script that compares a input GTF (-i) (typically "
                                     "reads in GTF form) to a gene model summary (ES "
                                     "GTF) (-es). Creates exon segment patterns (ESP) which are "
                                     "binary patterns indicating which of a gene's exon "
                                     "segments a transcript has exons within. Requires ER GTF "
                                     "(-er) to include ER delineation in the notation of the ESP."
                                     "Outputs two files: "
                                     "to desired output directory (-o). "
                                     "1. A list of transcripts and their ESPs. 2. A flag file "
                                     "that indicates which of the gene's exon segments the "
                                     "transcript has. Script will also output a list of genes"
                                     "that only appear in one GTF. There is also an option to add an "
                                     "output prefix (-x).")

    # INPUT
    parser.add_argument(
        "-i",
        "--input-gtf",
        dest="inFile",
        required=True,
        help="Location of input GTF"
    )

    parser.add_argument(
        "-es",
        "--es-gtf",
        dest="esFile",
        required=True,
        help="Location of ES GTF"
    )

    parser.add_argument(
        "-er",
        "--er-gtf",
        dest="erFile",
        required=True,
        help="Location of ER GTF"
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


def flagESPStructure(inDf):

    patternSeekDf = inDf.copy()

    # List of test patterns for dev
    # erpLst = [
    #     '1'*22,
    #     '1'*23,
    #     '0'*21+'1',
    #     '0'*22+'1',
    #     '0'*19+'1'*3,
    #     '0'*19+'1'*4,
    #     '0'*10 + '101' + '1'*9,
    #     '1' + '0'*21 + '1',
    #     '1110111111101111101100',
    #     '11101111111011111011001',
    #     '0000000000000000111111',
    #     '00000000000000001111111',
    #     '0000000000000000000001',
    #     '00000000000000000000011',
    #     '1111111110000000000000',
    #     '1000000000000000000000',
    #     '0000000011110000000000',
    #     '00000000111100000000001',
    #     '00000000100000000000001',
    #     '0000000010000000000000'
    # ]

    # geneLst = ['FBgn0004652'] * len(erpLst)
    # strandLst = ['-'] * len(erpLst)

    # patternSeekDf = pd.DataFrame({
    #     'geneID': geneLst,
    #     'ERP': erpLst,
    #     'strand': strandLst
    # })

    # Pattern discernment for describing ERPs!
    patternSeekDf['patternSeek'] = patternSeekDf['ESP'].str.split(
        '_').str[1]

    # Pattern discernment!
    # 1. flag transcripts with all exon segments in the gene and no reference exon segments
    patternSeekDf['flagNoSkip'] = patternSeekDf['patternSeek'].apply(
        lambda x: 1 if all(char == '1' for char in x) else 0)

    patternSeekDf['flagNovel'] = patternSeekDf['patternSeek'].apply(
        lambda x: 1 if all(char == '0' for char in x) else 0)

    # 2. flag transcripts with an exon skip (one missing ES between two present ESs)
    patternSeekDf['flagESSkip'] = patternSeekDf.apply(
        lambda x: 1 if re.search('(?<=1)+0+(?=1)+', x['patternSeek']) is not None else 0, axis=1)

    # 3. 5' and 3' fragment (compared to the gene)
    patternSeekDf['flag5pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            "^1+0+$", x['patternSeek']) is not None else 0, axis=1)

    patternSeekDf['flag3pFragment'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^0+1+$', x['patternSeek']) is not None else 0, axis=1)

    # 4. internal fragment
    patternSeekDf['flagIntrnlFrgmnt'] = patternSeekDf.apply(
        lambda x: 1 if re.search('^0+1+0+$', x['patternSeek']) is not None else 0, axis=1)

    # 5. first/last ES present
    patternSeekDf['flagFirstES'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^1', x['patternSeek']) is not None else 0, axis=1)

    patternSeekDf['flagLastES'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '1$', x['patternSeek']) is not None else 0, axis=1)

    return patternSeekDf


def main():

    # NOTE: This script has the **exact same logic** as the ERP version.
    # ## ^ NOT ANYMORE
    # esFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc_es.gtf"
    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_noMultiGene.gtf"
    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc_er.gtf"

    inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset.gtf"
    esFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset_es.gtf"
    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset_er.gtf"

    # outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # esFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_es.gtf"
    # # dataFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_noMultiGene.gtf"
    # dataFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    # erFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"

    # prefix = "test_prefix"
    sampleID = None
    prefix = None

    inFile = args.inFile
    esFile = args.esFile
    erFile = args.erFile
    outdir = args.outdir
    prefix = args.prefix
    sampleID = args.sampleID

    alphatic = time.perf_counter()

    # Read in both GTFs and subset them to genes that are in both GTF
    inGtfDf = trand.io.read_exon_data_from_file(inFile)
    inESDf = trand.io.read_exon_data_from_file(esFile)

    uniqGtfGeneSet = set(inGtfDf['gene_id'])
    uniqRefGeneSet = set(inESDf['gene_id'])

    # Store genes only in one GTF for later output
    inputOnlyGnLst = list(uniqGtfGeneSet - uniqRefGeneSet)
    refOnlyGnLst = list(uniqRefGeneSet - uniqGtfGeneSet)

    genesInBoth = list(uniqRefGeneSet.intersection(uniqGtfGeneSet))

    inDf = inGtfDf[inGtfDf['gene_id'].isin(genesInBoth)].copy()
    geneESDf = inESDf[inESDf['gene_id'].isin(genesInBoth)].copy()

    inERDf = trand.io.read_exon_data_from_file(erFile)
    geneERDf = inERDf[inERDf['gene_id'].isin(genesInBoth)].copy()

    # Clean up ES GTF dataframe (geneDf)
    geneESDf = geneESDf[['gene_id', 'seqname',
                         'start', 'end', 'strand']].copy()
    geneESDf = geneESDf.sort_values(
        ['seqname', 'gene_id', 'start'], ignore_index=True)

    geneERDf = geneERDf[['gene_id', 'seqname',
                         'start', 'end', 'strand']].copy()
    geneERDf = geneERDf.sort_values(
        ['seqname', 'gene_id', 'start'], ignore_index=True)

    # Check that each gene is only on one strand (don't know why they wouldn't be)
    singleStrandGene = geneESDf.groupby('gene_id').agg(
        set)['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandGene.all():
        print("There are genes belonging to more than one strand. Quitting.")
        quit()

    # Assign each exon in the ER GTF its ER ID
    geneERDf['ER'] = geneERDf['gene_id'] + ':ER' + \
        (geneERDf.groupby('gene_id').cumcount() + 1).astype(str)

    erDct = geneERDf.set_index('ER').to_dict(orient='index')
    geneERDct = dict(geneERDf.groupby('gene_id').apply(
        lambda x: sorted(set(x['ER']), key=lambda x: int(x.split("ER")[1]))))

    # Find overlapping exon region for every exon segment!
    rowDct = geneESDf.to_dict('records')

    for row in rowDct:

        gene = row['gene_id']
        esStart = row['start']
        esEnd = row['end']

        geneERLst = geneERDct[gene]

        for erID in geneERLst:

            erInfo = erDct[erID]

            erStart = erInfo['start']
            erEnd = erInfo['end']

            if esStart <= erEnd and erStart <= esEnd:
                row['ER'] = erID

    # Seems like it worked!
    geneERESDf = pd.DataFrame(rowDct)

    # Create ES IDs
    geneERESDf['ES'] = geneERESDf['ER'] + ':ES' + \
        (geneERESDf.groupby(['gene_id', 'ER']).cumcount() + 1).astype(str)

    if geneERESDf['ER'].isnull().any():
        raise Exception(
            "An error occurred when pairing exon segments to exon regions.")

    # Create a dictionary of genes/ERs and their ESs. Sort ESIDs to be in numerical order (matches 5'->3' relative to + strand)
    preDictDf = geneERESDf.groupby(['gene_id', 'ER']).apply(lambda x: sorted(
        set(x['ES']), key=lambda x: int(x.split("ES")[1]))).reset_index()

    preDictDf.columns = ['gene_id', 'ER', 'ES']

    geneERESDct = dict(dict())
    for row in preDictDf.to_dict('records'):

        gene = row['gene_id']
        ER = row['ER']
        ES = row['ES']

        if gene not in geneERESDct:
            geneERESDct[gene] = {ER: ES}

        geneERESDct[gene][ER] = ES

    geneESDct = dict(geneERESDf.groupby(['gene_id']).apply(lambda x: sorted(
        set(x['ES']), key=lambda x: int(x.split("ES")[1]))))

    # TODO: CHECK THAT ALL SETS ARE OF SIZE ONE
    # Create a dictionary of ESs and their information
    esDf = geneERESDf.groupby('ES').agg('first')
    esDf['length'] = esDf['end'] - esDf['start']
    esDct = esDf.to_dict(orient='index')

    inDf['numExon'] = inDf.groupby('transcript_id')[
        'transcript_id'].transform('count')

    inDf['dataOnlyExon'] = np.nan

    records = inDf.to_dict('records')

    for row in records:

        gene = row['gene_id']
        # jxnHash = row['transcript_id']

        matchingESIDLst = []

        # if gene in geneDct.keys():
        for esID in geneESDct.get(gene):
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
        raise Exception(
            "There are transcripts belonging to more than one strand. Quitting.")
    else:
        xscriptESDf['strand'] = xscriptESDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrXscript = xscriptESDf['seqname'].apply(lambda x: len(x) == 1)

    if not singleChrXscript.all():
        raise Exception(
            "There are transcripts belonging to more than one strand. Quitting.")
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
    esLst = []
    flagLst = []
    lngthLst = []

    # gene = "FBgn0287617"
    # transcript = "01b14131049c89b4265a67ceb9f6305e84bd050fd97b0efec4a51fcd3720e4bf"
    # strand = "+"
    # seqname = "2L"

    patternDct = dict()
    for gene, transcript, strand, seqname in loopLst:

        xscriptESSet = xscriptESDct.get(transcript)

        pttrnLst = []
        esIDLst = []
        # loop through ER then ES here. allows for separation of ERs in ESP notation
        for ER in list(geneERESDct.get(gene).keys()):

            for ES in geneERESDct.get(gene).get(ER):
                num = 1 if ES in xscriptESSet else 0
                pttrnLst.append(num)

                if ES in xscriptESSet:
                    esIDLst.append(ES)

            if ER != list(geneERESDct.get(gene).keys())[-1]:
                pttrnLst.append("-")
                # esIDLst.append("_")

        # pttrnLst = ["1" if ES in xscriptESSet else "0" for ES in geneESLst]
        # esIDLst = [ES for ES in geneESLst if ES in xscriptESSet]

        if strand == "-":
            pttrnLst.reverse()
            esIDLst.reverse()

        pattern = strand + "_" + ''.join(map(str, pttrnLst))
        patternESID = strand + "_" + "_".join(esIDLst)

        patternDct[transcript] = [pattern, patternESID, gene]

        geneESLst = geneESDct.get(gene)
        for exonSegment in geneESLst:

            if exonSegment in xscriptESSet:
                flag = 1
            else:
                flag = 0

            xscriptLst.append(transcript)
            geneLst.append(gene)
            esLst.append(exonSegment)
            flagLst.append(flag)
            lngthLst.append(esDct[exonSegment]['length'])

        if dataOnlyExonDct[transcript]:

            for exonRegion in dataOnlyExonDct[transcript]:
                xscriptLst.append(transcript)
                geneLst.append(gene)
                esLst.append(exonSegment)
                flagLst.append(1)

                # TODO: CHANGE THIS IF THE DATA ONLY EXON FORMAT CHANGES
                startNEnd = exonRegion.split(':')[1].split('_')
                exonLngth = int(startNEnd[1]) - int(startNEnd[0])

                lngthLst.append(exonLngth)

    outFlagDf = pd.DataFrame({
        'jxnHash': xscriptLst,
        'ES': esLst,
        'flagES': flagLst,
        'lengthES': lngthLst,
        'geneID': geneLst,
    })

    # Making pattern output file
    pttrnInfo = [(xscript, *info) for xscript, info in patternDct.items()]
    patternDf = pd.DataFrame(pttrnInfo, columns=[
        'transcript_id', 'ESP', 'patternES_ID', 'geneID'])

    outPatternDf = pd.merge(xscriptESDf, patternDf, on=[
        'transcript_id'], how='outer', indicator='merge_check')
    outPatternDf.rename(columns={'transcript_id': 'jxnHash'}, inplace=True)

    if not (outPatternDf['merge_check'] == 'both').all():
        raise Exception(
            "Something went wrong. Merge of patterns and xscript information failed.")

    outPatternDf['flagDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: len(x) != 0).astype(int)

    # TODO: no numES in ESP file
    outPatternDf['numDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(len)

    outPatternDf['dataOnlyES_ID'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    # Add ESP flags to output
    # outPatternDf = flagESPStructure(inDf=outPatternDf)

    # Output
    outFlagDf = outFlagDf.sort_values(by=['geneID', 'jxnHash'])
    outPatternDf = outPatternDf.sort_values(by=['geneID', 'jxnHash'])

    patternColLst = [
        'jxnHash',
        'geneID',
        'ESP',
        'numExon',
        'flagDataOnlyExon',
        'numDataOnlyExon',
        'dataOnlyES_ID',
        # 'flagNoSkip',
        # 'flagNovel',
        # 'flagERSkip',
        # 'flag5pFragment',
        # 'flag3pFragment',
        # 'flagIntrnlFrgmnt',
        # 'flagFirstER',
        # 'flagLastER',
        'patternES_ID'
    ]

    outPatternDf = outPatternDf[patternColLst]

    if sampleID:
        outFlagDf['sampleID'] = sampleID
        outPatternDf['sampleID'] = sampleID

    # Do not uncomment. Will probably crash the script.
    # wideDf = pd.pivot_table(outDf, values='flagES', index=['jxnHash','geneID'], columns='exonRegion', fill_value=0)

    esName = os.path.splitext(os.path.basename(esFile))[0]
    inName = os.path.splitext(os.path.basename(inFile))[0]

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    espFile = outPrefix + "{}_vs_{}_infoESP.csv".format(esName, inName)

    outPatternDf.to_csv(
        espFile, index=False)

    flagFile = outPrefix + "{}_vs_{}_flagES.csv".format(esName, inName)

    outFlagDf.to_csv(
        flagFile, index=False)

    if refOnlyGnLst:
        pd.Series(refOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_es_only_genes.txt".format(esName, inName), index=False, header=False)

    if inputOnlyGnLst:
        pd.Series(inputOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_input_only_genes.txt".format(esName, inName), index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

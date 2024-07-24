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
                                     "reads in GTF form) to an exon region (ER) "
                                     "GTF (-er). Creates exon region patterns (ERP) which are "
                                     "binary patterns indicating which of a gene's exon "
                                     "regions a transcript has exons within. Outputs two files: "
                                     "to desired output directory (-o). "
                                     "1. A list of transcripts and their ERPs. 2. A flag file "
                                     "that indicates which of the gene's exon regions the "
                                     "transcript has. Script will also output a list of genes"
                                     "that only appear one GTF. There is also an option to add an "
                                     "output prefix (-x). Note that there are two forms of "
                                     "these ERPs: one in which the exon regions "
                                     "match the way they would appear in the DNA (takes strand into account "
                                     "and one in which they match the order of the GTF. ")

    # INPUT
    parser.add_argument(
        "-er",
        "--er-gtf",
        dest="erFile",
        required=True,
        help="Location of ER GTF"
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

    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_data.gtf"
    # outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"

    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_2_dmel6_uniq_jxnHash_sexDet.gtf"

    # erFile = "//exasmb.rc.ufl.edu/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "//exasmb.rc.ufl.edu/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_data.gtf"
    # outdir = 'C://Users/knife/Desktop/Code Dumping Ground/mcintyre'

    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset_er.gtf"
    dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset.gtf"

    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc_sexDetSubset_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/gffcompare_read_aln_ujc/dyak_F_2_dyak2_ujc_updGeneID.gtf"

    outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # prefix = "test_prefix"
    sampleID = "testSampleID"
    prefix = None

    erFile = args.erFile
    dataFile = args.dataFile
    outdir = args.outdir
    prefix = args.prefix
    sampleID = args.sampleID

    alphatic = time.perf_counter()

    inGeneDf = trand.io.read_exon_data_from_file(erFile)
    inDataDf = trand.io.read_exon_data_from_file(dataFile)

    uniqDataGeneSet = set(inDataDf['gene_id'])
    uniqRefGeneSet = set(inGeneDf['gene_id'])

    refOnlyGnLst = list(uniqRefGeneSet - uniqDataGeneSet)
    dataOnlyGnLst = list(uniqDataGeneSet - uniqRefGeneSet)

    genesInBoth = list(uniqRefGeneSet.intersection(uniqDataGeneSet))

    geneDf = inGeneDf[inGeneDf['gene_id'].isin(genesInBoth)].copy()
    dataDf = inDataDf[inDataDf['gene_id'].isin(genesInBoth)].copy()

    geneDf = geneDf[['gene_id', 'seqname', 'start', 'end', 'strand']].copy()
    geneDf = geneDf.sort_values(
        ['seqname', 'gene_id', 'start'], ignore_index=True)

    singleStrandGene = geneDf.groupby('gene_id').agg(
        set)['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandGene.all():
        print("There are genes belonging to more than one strand. Quitting.")
        quit()

    geneDf['ER'] = geneDf['gene_id'] + ':ER' + \
        (geneDf.groupby('gene_id').cumcount() + 1).astype(str)

    geneDct = dict(geneDf.groupby('gene_id').apply(
        lambda x: sorted(set(x['ER']), key=lambda x: int(x.split("ER")[1]))))

    # TODO: CHECK THAT ALL SETS ARE OF SIZE ONE
    erDf = geneDf.groupby('ER').agg('first')
    erDf['length'] = erDf['end'] - erDf['start']
    erDct = erDf.to_dict(orient='index')

    # row = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+1])
    # row2 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+2])
    # row3 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+3])
    # yourBoat = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+4])
    # dataDf = pd.concat([dataDf,row,row2,row3,yourBoat])

    dataDf['numExon'] = dataDf.groupby('transcript_id')[
        'transcript_id'].transform('count')

    dataDf['dataOnlyExon'] = np.nan

    records = dataDf.to_dict('records')

    for row in records:

        gene = row['gene_id']
        # jxnHash = row['transcript_id']

        matchingERIDLst = []

        if gene in geneDct.keys():
            for erID in geneDct.get(gene):
                # print(erID)
                erInfo = erDct.get(erID)
                # print(erInfo)

                # print("looping...")

                if max(row['start'], erInfo['start']) < min(row['end'], erInfo['end']):
                    # print(row)
                    # print(erID)
                    # print(erInfo)

                    matchingERIDLst.append(erID)

            if matchingERIDLst:
                row['ER'] = matchingERIDLst
            else:
                row['dataOnlyExon'] = "{}:{}_{}".format(
                    gene, row['start'], row['end'])

    dataWithERDf = pd.DataFrame(records)

    # TODO: May not be true but makes sense right? if there is an exon that overlaps two exon regions
    # it may not be true biological IR but its definitely overlapping an intron...
    dataWithERDf['flagIR'] = dataWithERDf['ER'].apply(
        lambda x: x if not type(x) is list else 1 if len(x) > 1 else 0)

    dataWithERDf['IRER'] = dataWithERDf.apply(
        lambda x: tuple(x['ER']) if x['flagIR'] == 1 else np.nan, axis=1)

    dataWithERDf['numIREvent'] = dataWithERDf.groupby(
        'transcript_id')['flagIR'].transform('sum')

    intmdDf = dataWithERDf[['gene_id', 'transcript_id', 'ER',
                            'dataOnlyExon', 'flagIR', 'numIREvent', 'IRER', 'numExon', 'strand']]
    intmdDf = intmdDf.explode('ER')

    xscriptERDf = intmdDf.groupby('transcript_id').agg({
        'ER': lambda x: set(x.dropna()),
        'numExon': max,
        'flagIR': max,
        'dataOnlyExon': lambda x: set(x.dropna()),
        'IRER': lambda x: set(tuple(sum(x.dropna(), ()))),
        'strand': set,
        'numIREvent': max
    }).reset_index()

    # Check for no multi-strand transcripts
    singleStrandXscript = xscriptERDf['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandXscript.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        xscriptERDf['strand'] = xscriptERDf['strand'].apply(
            lambda x: list(x)[0])

    # Accounts for situations where transcripts have multiple IR events

    xscriptERDct = dict(zip(xscriptERDf['transcript_id'], xscriptERDf['ER']))

    dataOnlyExonDct = dict(
        zip(xscriptERDf['transcript_id'], xscriptERDf['dataOnlyExon']))

    # dataWithERDf = dataWithERDf[
    #     (dataWithERDf['transcript_id'] ==
    #       "34a6dd0389208ff91783b8ec557da2427785a0988ffa8243635e75c13b7f334c")
    #     | (dataWithERDf['transcript_id'] == "068509f23790d261052860383c257e94c8d917c9c67a0140875242cd7302b738")]

    loopLst = [tuple(x) for x in dataWithERDf[[
        'gene_id', 'transcript_id', 'strand']].drop_duplicates().to_records(index=False)]

    xscriptLst = []
    geneLst = []
    erLst = []
    flagLst = []
    strandLst = []
    lngthLst = []

    patternDct = dict()
    for gene, transcript, strand in loopLst:

        geneERLst = geneDct.get(gene)
        xscriptERSet = xscriptERDct.get(transcript)

        pttrnLst = ["1" if ER in xscriptERSet else "0" for ER in geneERLst]
        erIDLst = [ER for ER in geneERLst if ER in xscriptERSet]

        if strand == "-":
            pttrnLst.reverse()
            erIDLst.reverse()

        pattern = strand + "_" + ''.join(map(str, pttrnLst))
        patternERID = strand + "_" + "_".join(erIDLst)

        patternDct[transcript] = [pattern, patternERID, gene]

        for exonRegion in geneERLst:

            if exonRegion in xscriptERSet:
                flag = 1
            else:
                flag = 0

            xscriptLst.append(transcript)
            geneLst.append(gene)
            erLst.append(exonRegion)
            flagLst.append(flag)
            strandLst.append(strand)
            lngthLst.append(erDct[exonRegion]['length'])

        if dataOnlyExonDct[transcript]:

            for exonRegion in dataOnlyExonDct[transcript]:
                xscriptLst.append(transcript)
                geneLst.append(gene)
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
        'strand': strandLst,
        'ER': erLst,
        'flagER': flagLst,
        'lengthER': lngthLst
    })

    # Making pattern output file
    pttrnInfo = [(xscript, *info) for xscript, info in patternDct.items()]
    patternDf = pd.DataFrame(pttrnInfo, columns=[
        'transcript_id', 'ERP', 'patternERID', 'geneID'])

    outPatternDf = pd.merge(xscriptERDf, patternDf, on=[
        'transcript_id'], how='outer', indicator='merge_check')
    outPatternDf.rename(columns={'transcript_id': 'jxnHash'}, inplace=True)

    if not (outPatternDf['merge_check'] == 'both').all():
        raise Exception(
            "Something went wrong. Merge of patterns and xscript information failed.")

    outPatternDf['flagDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: len(x) != 0)

    # TODO: no numER in ERP file
    outPatternDf['numER'] = outPatternDf['ER'].apply(len)
    outPatternDf['numDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(len)

    # TODO: rename reverseIR
    outPatternDf['flagReverseIR'] = outPatternDf.apply(
        lambda x: 1 if x['numExon'] > x['numER'] + x['numDataOnlyExon'] else 0, axis=1)

    outPatternDf['IRER'] = outPatternDf['IRER'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    outPatternDf['dataOnlyERID'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    outFlagDf = outFlagDf.sort_values(by=['geneID', 'jxnHash'])
    outPatternDf = outPatternDf.sort_values(by=['geneID', 'jxnHash'])

    outPatternDf = outPatternDf[['jxnHash', 'geneID', 'strand', 'ERP', 'patternERID', 'numExon',
                                 'numDataOnlyExon', 'dataOnlyERID', 'flagIR', 'numIREvent', 'IRER',
                                 'flagReverseIR']]

    if sampleID:
        outFlagDf['sampleID'] = sampleID
        outPatternDf['sampleID'] = sampleID

    # Do not uncomment. Will probably crash the script.
    # wideDf = pd.pivot_table(outDf, values='flag_ER', index=['jxnHash','geneID'], columns='exonRegion', fill_value=0)

    # Output
    erName = os.path.splitext(os.path.basename(erFile))[0]
    dataName = os.path.splitext(os.path.basename(dataFile))[0]

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    erpFile = outPrefix + "{}_vs_{}_ERP.csv".format(erName, dataName)

    outPatternDf.to_csv(
        erpFile, index=False)

    flagFile = outPrefix + "{}_vs_{}_flagER.csv".format(erName, dataName)

    outFlagDf.to_csv(
        flagFile, index=False)

    if refOnlyGnLst:
        pd.Series(refOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_anno_only_genes.txt".format(erName, dataName), index=False, header=False)

    if dataOnlyGnLst:
        pd.Series(dataOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_data_only_genes.txt".format(erName, dataName), index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

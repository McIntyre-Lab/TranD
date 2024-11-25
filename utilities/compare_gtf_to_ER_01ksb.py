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
    parser = argparse.ArgumentParser(description="A script that compares an input GTF (-i) (typically "
                                     "reads in GTF form) to a gene model summary (ER "
                                     "GTF) (-er). Creates exon region patterns (ERP) which are "
                                     "binary patterns indicating which of a gene's exon "
                                     "regions a read has exons within. Outputs two files: "
                                     "to desired output directory (-o). "
                                     "1. A list of transcripts and their ERPs. 2. A flag file "
                                     "that indicates which of the gene's exon regions the "
                                     "transcript has. Script will also output a list of genes"
                                     "that only appear one GTF. There is also an option to add an "
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


def flagERPStructure(inDf):

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
    patternSeekDf['patternSeek'] = patternSeekDf['ERP'].str.split(
        '_').str[1]

    # 1. flag transcripts with all exon regions in the gene and no reference exon regions
    patternSeekDf['flagNoSkip'] = patternSeekDf['patternSeek'].apply(
        lambda x: 1 if all(char == '1' for char in x) else 0)

    patternSeekDf['flagNovel'] = patternSeekDf['patternSeek'].apply(
        lambda x: 1 if all(char == '0' for char in x) else 0)

    # 2. flag transcripts with an exon skip (one missing ER between two present ERs)
    patternSeekDf['flagERSkip'] = patternSeekDf.apply(
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

    # 5. first/last ER present
    patternSeekDf['flagFirstER'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '^1', x['patternSeek']) is not None else 0, axis=1)

    patternSeekDf['flagLastER'] = patternSeekDf.apply(
        lambda x: 1 if re.search(
            '1$', x['patternSeek']) is not None else 0, axis=1)

    return patternSeekDf


def main():

    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_data.gtf"
    # # outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"

    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/mel_2_dmel6_uniq_jxnHash_sexDet.gtf"

    # erFile = "//exasmb.rc.ufl.edu/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "//exasmb.rc.ufl.edu/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_data.gtf"
    # outdir = 'C://Users/knife/Desktop/Code Dumping Ground/mcintyre'

    inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset.gtf"
    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_sexDetSubset_er.gtf"

    # inFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/dyak_data_2_dyak2_ujc_noMultiGene.gtf"
    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc_er.gtf"

    outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # prefix = "test_prefix"
    sampleID = "testSampleID"
    prefix = None

    erFile = args.erFile
    inFile = args.inFile
    outdir = args.outdir
    prefix = args.prefix
    sampleID = args.sampleID

    alphatic = time.perf_counter()

    # Read in both GTFs and subset them to genes that are in both GTF
    inGtfDf = trand.io.read_exon_data_from_file(inFile)
    inERDf = trand.io.read_exon_data_from_file(erFile)

    uniqGtfGeneSet = set(inGtfDf['gene_id'])
    uniqRefGeneSet = set(inERDf['gene_id'])

    # Store genes only in one GTF for later output
    inputOnlyGnLst = list(uniqGtfGeneSet - uniqRefGeneSet)
    refOnlyGnLst = list(uniqRefGeneSet - uniqGtfGeneSet)

    genesInBoth = list(uniqRefGeneSet.intersection(uniqGtfGeneSet))

    inDf = inGtfDf[inGtfDf['gene_id'].isin(genesInBoth)].copy()
    geneERDf = inERDf[inERDf['gene_id'].isin(genesInBoth)].copy()

    # geneDf = inGeneDf[inGeneDf['gene_id'].isin(["LOC120456871"])].copy()
    # dataDf = inDataDf[inDataDf['gene_id'].isin(["LOC120456871"])].copy()

    # Clean up ER GTF dataframe (geneDf)
    geneERDf = geneERDf[['gene_id', 'seqname',
                         'start', 'end', 'strand']].copy()
    geneERDf = geneERDf.sort_values(
        ['seqname', 'gene_id', 'start'], ignore_index=True)

    # Check that each gene is only on one strand (don't know why they wouldn't be)
    singleStrandGene = geneERDf.groupby('gene_id').agg(
        set)['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandGene.all():
        print("There are genes belonging to more than one strand. Quitting.")
        quit()

    # Assign each exon in the ER GTF its ER ID
    geneERDf['ER'] = geneERDf['gene_id'] + ':ER' + \
        (geneERDf.groupby('gene_id').cumcount() + 1).astype(str)

    # Create a dictionary of genes and their ERs. Sort ERIDs to be in numerical order (matches 5'->3' relative to + strand)
    geneDct = dict(geneERDf.groupby('gene_id').apply(
        lambda x: sorted(set(x['ER']), key=lambda x: int(x.split("ER")[1]))))

    # TODO: CHECK THAT ALL SETS ARE OF SIZE ONE
    # Create a dictionary of ERs and their information
    erDf = geneERDf.groupby('ER').agg('first')
    erDf['length'] = erDf['end'] - erDf['start']
    erDct = erDf.to_dict(orient='index')

    # row = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+1])
    # row2 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+2])
    # row3 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+3])
    # yourBoat = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+4])
    # dataDf = pd.concat([dataDf,row,row2,row3,yourBoat])

    # Prepare data GTF/Dataframe for assigning ER
    inDf['numExon'] = inDf.groupby('transcript_id')[
        'transcript_id'].transform('count')

    inDf['dataOnlyExon'] = np.nan

    records = inDf.to_dict('records')

    # Loop through every exon in the data and create a list of the ER(s) it overlaps
    for row in records:

        gene = row['gene_id']
        # jxnHash = row['transcript_id']

        matchingERIDLst = []

        # if gene in geneDct.keys():
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

        # If the exon does not overlap any annotated ERs, add it to a list of "dataOnlyExon"
        if matchingERIDLst:
            row['ER'] = matchingERIDLst
        else:
            row['dataOnlyExon'] = "{}:{}_{}".format(
                gene, row['start'], row['end'])

    # Create new dataframe where every exon is assigned to its ER(s)
    dataWithERDf = pd.DataFrame(records)

    # flag exons with IR (one exon overlaps multiple ERs)
    dataWithERDf['flagIR'] = dataWithERDf['ER'].apply(
        lambda x: x if not type(x) is list else 1 if len(x) > 1 else 0)
    dataWithERDf['flagIR'] = dataWithERDf['flagIR'].fillna(0).astype(int)

    # List the ERs involved in the IR event
    dataWithERDf['IR_ER'] = dataWithERDf.apply(
        lambda x: tuple(x['ER']) if x['flagIR'] == 1 else np.nan, axis=1)

    # Give total number of times there is IR in that transcript (number of exons that have IR)
    dataWithERDf['numIREvent'] = dataWithERDf.groupby(
        'transcript_id')['flagIR'].transform('sum')

    # Create an intermediate DF based on new data exon DF
    intmdDf = dataWithERDf[['seqname', 'gene_id', 'transcript_id', 'ER',
                            'dataOnlyExon', 'flagIR', 'numIREvent', 'IR_ER', 'numExon', 'strand']]
    intmdDf = intmdDf.explode('ER')

    # Create a DF unique on transcript that lists all ERs associated with transcript
    xscriptERDf = intmdDf.groupby('transcript_id').agg({
        'ER': lambda x: set(x.dropna()),
        'numExon': max,
        'flagIR': max,
        'dataOnlyExon': lambda x: set(x.dropna()),
        'IR_ER': lambda x: set(tuple(sum(x.dropna(), ()))),
        'strand': set,
        'numIREvent': max,
        'seqname': set
    }).reset_index()

    # Check for no multi-strand or multi-chr transcripts post-grouping
    singleStrandXscript = xscriptERDf['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandXscript.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        xscriptERDf['strand'] = xscriptERDf['strand'].apply(
            lambda x: list(x)[0])

    singleChrXscript = xscriptERDf['seqname'].apply(lambda x: len(x) == 1)

    if not singleChrXscript.all():
        print("There are transcripts belonging to more than one strand. Quitting.")
        quit()
    else:
        xscriptERDf['seqname'] = xscriptERDf['seqname'].apply(
            lambda x: list(x)[0])

    # Create a dict that lists every annotated ER associated with a transcript
    # Accounts for situations where transcripts have multiple IR events
    xscriptERDct = dict(zip(xscriptERDf['transcript_id'], xscriptERDf['ER']))

    # Create a dict that lists every non-annotated exon assocaited with a transcript
    dataOnlyExonDct = dict(
        zip(xscriptERDf['transcript_id'], xscriptERDf['dataOnlyExon']))

    # dataWithERDf = dataWithERDf[
    #     (dataWithERDf['transcript_id'] ==
    #       "34a6dd0389208ff91783b8ec557da2427785a0988ffa8243635e75c13b7f334c")
    #     | (dataWithERDf['transcript_id'] == "068509f23790d261052860383c257e94c8d917c9c67a0140875242cd7302b738")]

    # Create a list of transcripts to loop through using the data GTF
    loopLst = [tuple(x) for x in dataWithERDf[[
        'gene_id', 'transcript_id', 'strand', 'seqname']].drop_duplicates().to_records(index=False)]

    # Loop through every data transcript.
    # Compare the ERs in the transcript to every ER in the gene and flag present ERs.
    # Importantly, the ERs are looped through the gene's list of ERs, which are in order from
    # 5' -> 3' due to the previous sort.
    # Create the ERP

    xscriptLst = []
    geneLst = []
    erLst = []
    flagLst = []
    lngthLst = []

    patternDct = dict()
    for gene, transcript, strand, seqname in loopLst:

        geneERLst = geneDct.get(gene)
        xscriptERSet = xscriptERDct.get(transcript)

        pttrnLst = ["1" if ER in xscriptERSet else "0" for ER in geneERLst]
        erIDLst = [ER for ER in geneERLst if ER in xscriptERSet]

        # Assure that ERPs are in transcription order in the pattern
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
            lngthLst.append(erDct[exonRegion]['length'])

        # Add any data only exons if present
        if dataOnlyExonDct[transcript]:

            for exonRegion in dataOnlyExonDct[transcript]:
                xscriptLst.append(transcript)
                geneLst.append(gene)
                erLst.append(exonRegion)
                flagLst.append(1)

                # TODO: CHANGE THIS IF THE DATA ONLY EXON FORMAT CHANGES
                startNEnd = exonRegion.split(':')[1].split('_')
                exonLngth = int(startNEnd[1]) - int(startNEnd[0])

                lngthLst.append(exonLngth)

    # Create flagER file using lists created above
    outFlagDf = pd.DataFrame({
        'jxnHash': xscriptLst,
        'ER': erLst,
        'flagER': flagLst,
        'lengthER': lngthLst,
        'geneID': geneLst,
    })

    # Making pattern output file
    pttrnInfo = [(xscript, *info) for xscript, info in patternDct.items()]
    patternDf = pd.DataFrame(pttrnInfo, columns=[
        'transcript_id', 'ERP', 'patternER_ID', 'geneID'])

    outPatternDf = pd.merge(xscriptERDf, patternDf, on=[
        'transcript_id'], how='outer', indicator='merge_check')
    outPatternDf.rename(columns={'transcript_id': 'jxnHash'}, inplace=True)

    if not (outPatternDf['merge_check'] == 'both').all():
        raise Exception(
            "Something went wrong. Merge of patterns and xscript information failed.")

    outPatternDf['flagDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: len(x) != 0).astype(int)

    outPatternDf['numDataOnlyExon'] = outPatternDf['dataOnlyExon'].apply(len)

    # TODO: rename reverseIR
    # reverseIR = when multiple exons fit into one exon region
    outPatternDf['flagReverseIR'] = outPatternDf.apply(
        lambda x: 1 if x['numExon'] > len(x['ER']) + x['numDataOnlyExon'] else 0, axis=1)

    outPatternDf['IR_ER'] = outPatternDf['IR_ER'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    outPatternDf['dataOnlyER_ID'] = outPatternDf['dataOnlyExon'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    outPatternDf['ERP_plus'] = outPatternDf['ERP'].astype(
        str) + '|' + outPatternDf['flagDataOnlyExon'].astype(str) + '|' + outPatternDf['flagIR'].astype(str)

    # Add ERP flags to output
    outPatternDf = flagERPStructure(inDf=outPatternDf)

    # Output
    outFlagDf = outFlagDf.sort_values(by=['geneID', 'jxnHash'])
    outPatternDf = outPatternDf.sort_values(by=['geneID', 'jxnHash'])

    patternColLst = [
        'jxnHash',
        'geneID',
        'ERP',
        'ERP_plus',
        'numExon',
        'flagDataOnlyExon',
        'numDataOnlyExon',
        'dataOnlyER_ID',
        'flagIR',
        'numIREvent',
        'IR_ER',
        'flagReverseIR',
        'flagNoSkip',
        'flagNovel',
        'flagERSkip',
        'flag5pFragment',
        'flag3pFragment',
        'flagIntrnlFrgmnt',
        'flagFirstER',
        'flagLastER',
        'patternER_ID'
    ]

    outPatternDf = outPatternDf[patternColLst]

    if sampleID:
        outFlagDf['sampleID'] = sampleID
        outPatternDf['sampleID'] = sampleID

    # Do not uncomment. Will probably crash the script.
    # wideDf = pd.pivot_table(outDf, values='flag_ER', index=['jxnHash','geneID'], columns='exonRegion', fill_value=0)

    erName = os.path.splitext(os.path.basename(erFile))[0]
    inName = os.path.splitext(os.path.basename(inFile))[0]

    if prefix:
        outPrefix = "{}/{}_".format(outdir, prefix)
    else:
        outPrefix = "{}/".format(outdir)

    erpFile = outPrefix + "{}_vs_{}_infoERP.csv".format(erName, inName)

    outPatternDf.to_csv(
        erpFile, index=False)

    flagFile = outPrefix + "{}_vs_{}_flagER.csv".format(erName, inName)

    outFlagDf.to_csv(
        flagFile, index=False)

    if refOnlyGnLst:
        pd.Series(refOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_er_only_genes.txt".format(erName, inName), index=False, header=False)

    if inputOnlyGnLst:
        pd.Series(inputOnlyGnLst).to_csv(
            outPrefix + "list_{}_vs_{}_input_only_genes.txt".format(erName, inName), index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

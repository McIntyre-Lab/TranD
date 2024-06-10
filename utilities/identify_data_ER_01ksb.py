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


def getOptions():
    """

    Function to store user input via argparse

    Returns
    -------
    args : ARGPARSE ARGUMENTS
            User input via argparse.

    """

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Placeholder")

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
    parser.add_argument("-o",
                        "--outdir",
                        dest="outdir",
                        required=True,
                        help="Output file path")

    parser.add_argument("-x",
                        "--prefix",
                        dest="prefix",
                        required=True,
                        help="Output prefix. Required.")

    args = parser.parse_args()
    return args


def main():

    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_data.gtf"
    outdir = "/nfshome/k.bankole/Desktop/test_folder"

    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/dmel_sexdet_fiveSpecies_vs_data/mel_2_dmel6_uniq_jxnHash_sexDet.gtf"

    erFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    dataFile = "//exasmb.rc.ufl.edu/blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_data.gtf"
    outdir = 'C://Users/knife/Desktop/Code Dumping Ground/mcintyre'

    prefix = "test_prefix"

    erFile = args.erFile
    dataFile = args.dataFile
    outdir = args.outdir
    prefix = args.prefix

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
    geneDf = geneDf.sort_values(['seqname', 'gene_id', 'start'])

    geneDf['ER'] = geneDf['gene_id'] + ':ER' + \
        (geneDf.groupby('gene_id').cumcount() + 1).astype(str)
    
    singleStrandGene = geneDf.groupby('gene_id').agg(set)['strand'].apply(lambda x: len(x) == 1)

    if not singleStrandGene.all():
        print("There are genes belonging to more than one strand. Quitting.")
        quit()
    
    # sorted(set(
    #     x['ER']), key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("exon_")[1])))
    
        
    geneDct = dict(geneDf.groupby('gene_id').apply(lambda x: 
        sorted(set(x['ER']), key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("exon_")[1])) 
        if (x['strand'] == "+").all() 
        else 
        sorted(set(x['ER']), key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("exon_")[1]),reverse=True)))
        
    # TODO: CHECK THAT ALL SETS ARE OF SIZE ONE
    erDct = geneDf.groupby('ER').agg('first').to_dict(orient='index')

    # row = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+1])
    # row2 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+2])
    # row3 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+3])
    # yourBoat = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+4])

    # dataDf = pd.concat([dataDf,row,row2,row3,yourBoat])

    dataDf['numExon'] = dataDf.groupby('transcript_id')[
        'transcript_id'].transform('count')

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
                    # print ("Yippee!")
                    # print(row)
                    # print(erID)
                    # print(erInfo)
                    matchingERIDLst.append(erID)

            if matchingERIDLst:
                row['ER'] = matchingERIDLst
            else:
                row['dataOnlyER'] = "{}:{}_{}".format(
                    gene, row['start'], row['end'])

    dataWithERDf = pd.DataFrame(records)

    # TODO: May not be true but makes sense right? if there is an exon that overlaps two exon regions
    # it may not be true biological IR but its definitely overlapping an intron...
    dataWithERDf['flagIR'] = dataWithERDf['ER'].apply(
        lambda x: x if not type(x) is list else 1 if len(x) > 1 else 0)
    dataWithERDf['IRERs'] = dataWithERDf.apply(
        lambda x: tuple(x['ER']) if x['flagIR'] == 1 else np.nan, axis=1)

    intmdDf = dataWithERDf[['gene_id', 'transcript_id', 'ER',
                            'dataOnlyER', 'flagIR', 'IRERs', 'numExon', 'strand']]
    intmdDf = intmdDf.explode('ER')

    xscriptERDf = intmdDf.groupby('transcript_id').agg({
        'ER': lambda x: set(x.dropna()),
        'numExon': max,
        'flagIR': max,
        'dataOnlyER': lambda x: set(x.dropna()),
        'IRERs': lambda x: set(tuple(sum(x.dropna(), ()))),
        'strand': set
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
    xscriptERDf['numIREvent'] = xscriptERDf.groupby(
        'transcript_id')['flagIR'].transform('sum')

    xscriptERDct = dict(zip(xscriptERDf['transcript_id'], xscriptERDf['ER']))

    dataOnlyERDct = dict(
        zip(xscriptERDf['transcript_id'], xscriptERDf['dataOnlyER']))

    loopLst = [tuple(x) for x in dataWithERDf[[
        'gene_id', 'transcript_id', 'strand']].drop_duplicates().to_records(index=False)]

    xscriptLst = []
    geneLst = []
    erLst = []
    flagLst = []
    strandLst = []

    binaryDct = dict()
    for gene, transcript, strand in loopLst:

        # gene = row['geneID']
        # transcript = row['jxnHash']

        geneERLst = geneDct.get(gene)
        xscriptERSet = xscriptERDct.get(transcript)

        binary = [1 if ER in xscriptERSet else 0 for ER in geneERLst]
        binary = ''.join(map(str, binary))
        binaryDct[transcript] = [binary, gene]

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

        if dataOnlyERDct[transcript]:

            for exonRegion in dataOnlyERDct[transcript]:
                xscriptLst.append(transcript)
                geneLst.append(gene)
                erLst.append(exonRegion)
                flagLst.append(1)
                strandLst.append(strand)


    outFlagDf = pd.DataFrame({
        'jxnHash': xscriptLst,
        'geneID': geneLst,
        'strand':strandLst,
        'exonRegion': erLst,
        'flag_ERPresent': flagLst
    })
    

    # Making pattern output file
    binaryInfo = [(xscript, *info) for xscript, info in binaryDct.items()]
    binaryDf = pd.DataFrame(binaryInfo, columns=[
                            'transcript_id', 'ERP', 'geneID'])


    outPatternDf = pd.merge(xscriptERDf, binaryDf, on=[
        'transcript_id'], how='outer', indicator='merge_check')
    outPatternDf.rename(columns={'transcript_id':'jxnHash'}, inplace=True)
    
    if not (outPatternDf['merge_check'] == 'both').all():
        print("Something went wrong")
        quit()
    
    outPatternDf['flagDataOnlyExon'] = outPatternDf['dataOnlyER'].apply(lambda x: len(x) != 0)
    
    outPatternDf['numER'] = outPatternDf['ER'].apply(len)
    outPatternDf['numDataOnlyER'] = outPatternDf['dataOnlyER'].apply(len)

    
    outPatternDf = outPatternDf[['jxnHash', 'geneID', 'strand', 'ERP', 'numExon',
                                 'numER', 'numDataOnlyER', 'flagIR', 'numIREvent', 'IRERs']]

    outPatternDf['flagReverseIR'] = outPatternDf.apply(
        lambda x: 1 if x['numExon'] > x['numER'] + x['numDataOnlyER'] else 0, axis=1)

    outPatternDf['IRERs'] = outPatternDf['IRERs'].apply(
        lambda x: '|'.join(x) if x else np.nan)

    outFlagDf = outFlagDf.sort_values(by=['geneID', 'jxnHash'])
    outPatternDf = outPatternDf.sort_values(by=['geneID', 'jxnHash'])

    # Do not uncomment. Will probably crash the script.
    # wideDf = pd.pivot_table(outDf, values='flag_ERPresent', index=['jxnHash','geneID'], columns='exonRegion', fill_value=0)

    # Output
    outPatternDf.to_csv(
        "{}/{}_er_vs_data_pattern_file.csv".format(outdir, prefix), index=False)

    outFlagDf.to_csv(
        "{}/{}_er_vs_data_flag_file.csv".format(outdir, prefix), index=False)

    if refOnlyGnLst:
        pd.Series(refOnlyGnLst).to_csv(
            "{}/list_{}_er_vs_data_reference_only_genes.txt".format(outdir, prefix), index=False, header=False)

    if dataOnlyGnLst:
        pd.Series(dataOnlyGnLst).to_csv(
            "{}/list_{}_er_vs_data_data_only_genes.txt".format(outdir, prefix), index=False, header=False)

    omegatoc = time.perf_counter()

    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

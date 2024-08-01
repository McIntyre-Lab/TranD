#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on Tue Sep 26 13:34:55 2023

@author: k.bankole
"""

"""

Identify the unique junction chains (UJC) of a GTF file and combine
transcripts into UJCs.

Created from a previous utility in TranD named consolidation

TranD version of the utility is referred to as 1.xx in versioning.

Version 2.3: Output file prefix is now just input file name; junction string now only contains an 
                underscore connecting all junction positions.
"""

# import pickle




import argparse
import time
import pandas as pd
import os
import csv
import hashlib
import copy
import sys
import os
import trand.io
def getOptions():
    """

    Function to store user input via argparse

    Returns
    -------
    args : ARGPARSE ARGUMENTS
            User input via argparse.

    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Identifies unique junction chains (UJCs) found within a "
                                     "GTF file and groups transcripts based on these UJCs. Outputs 3 files: "
                                     "1. An id file that is unique on jxnHash and provides the following "
                                     "information: gene, chr, strand, xscript, numJxn, start, end, jxnString "
                                     "and jxnHash. jxnStrings are in the format of chr_strand_start:end_start:end_start:end "
                                     "with monoexons as chr_strand_monoexon_start_end (with the start/end being for the exon)."
                                     "jxnHash is the unique identifier for the group of transcripts. 2. A summary"
                                     "file unique on transcript that more directly links a read/transcript to its "
                                     "jxnHash and jxnString. 3. A GTF file with representative transcript models "
                                     "for each group. The jxnHash will be the \'transcript_id\' for the group. "
                                     "Input: a GTF file (--gtf), and an output directory (--outdir). Output files will begin with "
                                     "the name of the input file. Allows the option to skip the output of the GTF file "
                                     "with representative transcript models. (--skip-gtf). Allows the option to output "
                                     "another key file with the number of transcripts per jxnHash counted (--count-ujc).")

    # INPUT
    parser.add_argument(
        "-g",
        "--gtf",
        dest="inGTF",
        required=True,
        help="Input a GTF file."
    )

    parser.add_argument(
        "-s",
        "--skip-gtf",
        dest="includeGTF",
        action="store_false",
        help="Use this argument to remove the output of a GTF with "
        "representative transcript models for each UJC."
        "Defaults to outputting the GTF."
    )

    parser.add_argument(
        "-c",
        "--count-ujc",
        dest="includeCnt",
        action="store_true",
        help="Use this argument to output a key file that counts"
        "the number of transcripts per UJC. Defaults to no output."
    )

    parser.add_argument(
        "-kg",
        "--keep-geneID",
        dest="keepGene",
        action="store_true",
        help="Use this argument to keep the geneID from the "
        "input transcript and apply it to the jxnHash. WARNING: THIS "
        "WILL ONLY WORK PROPERLY IF THERE IS ONLY ONE GENE IN THE GTF "
        "(since it is possible for one jxnHash to belong to more "
        "than one gene)."
    )

    parser.add_argument(
        "-ts",
        "--track-source",
        dest="trackSrc",
        action="store_true",
        help="Use this argument to keep the source from the "
        "input GTF. Useful if inputting a GTF with multiple catted samples. "
        "make sure the sample is in the source column of the GTF. Will create "
        "a sampleID column in xscript_link and ujc_count files. This option "
        "will fail if there are duplicate transcriptIDs with two different "
        "samples/sources."
    )

    # OUTPUT
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Location of output directory, must already exist."
    )

    # parser.add_argument(
    #         "-p",
    #         "--prefix",
    #         dest="prefix",
    #         required=True,
    #         help="Required prefix for the output file(s). Example: prefix_UJC_ID.csv"
    # )

    args = parser.parse_args()
    return args


def checkStrandAndChromosome(exonData):
    """

    Checks each strand and chromosome to see if there are genes with transcripts/
    exons on both strands/different chromosomes and removes them.

    Parameters
    ----------
    exonData : DATAFRAME
            A GTF converted to a DataFrame with exon data.

    Returns
    -------
    exonData : DATAFRAME
            The same input with genes removed if necessary.

    """

    geneGrps = exonData.groupby("gene_id")
    strandCheck = geneGrps["strand"].nunique()

    if (strandCheck > 1).any():
        badStrand = list(strandCheck[strandCheck > 1].index)
        for gene in badStrand:
            print("!!! WARNING: gene {} contains transcripts/exons on both strands - "
                  "removing.".format(gene))

        exonData = exonData[~exonData["gene_id"].isin(badStrand)]

    chrCheck = geneGrps["seqname"].nunique()
    if (chrCheck > 1).any():
        badChr = list(chrCheck[chrCheck > 1].index)
        for gene in badChr:
            print("!!! WARNING: gene {} contains transcripts/exons on difference chromosomes - "
                  "removing.".format(gene))

        exonData = exonData[~exonData["gene_id"].isin(badChr)]

    return exonData


def extractJunction(exonData, trackSrc=False):
    """

    Takes the exon data and extracts the locations of the junctions for each
    transcript. Outputs information on each transcript to a dictionary with
    the transcript as a key and a list of information in the following format:
            [[exonLst], transcript_id, gene_id, seqname, start, end, strand]


    Parameters
    ----------
    exonData : DATAFRAME
            A GTF converted to a DataFrame with exon data.

    Returns
    -------
    ujcDct : DICTIONARY {Transcript_id: [info]}
            A dictionary of transcripts keyed to their info.

    """

    exonDf = checkStrandAndChromosome(exonData=exonData)

    print("Number of transcripts: ", end="")
    print(len(exonDf['transcript_id'].unique()))

    print("Number of genes: ", end="")
    print(len(exonDf['gene_id'].unique()))

    # First, instead of grouping, then sorting
    # Sort by transcript -> sort by start. the whole dataframe
    sortedDf = exonDf.sort_values(
        by=['transcript_id', 'start']).reset_index(drop=True)

    ujcDct = {}
    # Value Legend:
    # 0 = exons
    # 1 = xscript
    # 2 = gene
    # 3 = seqname
    # 4 = start
    # 5 = end
    # 6 = strand
    # 7 = source
    # 8 = junction string
    # 9 = numJxn

    for row in sortedDf.to_dict('records'):
        xscript = row['transcript_id']

        seqname = row['seqname']
        strand = row['strand']
        geneID = row['gene_id']

        start = row['start']
        end = row['end']

        if trackSrc:
            source = row['source']
        else:
            source = ""

        if xscript in ujcDct.keys():
            info = ujcDct[xscript]
            exonLst = info[0]
            oldStart = info[4]
            oldEnd = info[5]

            if oldStart > start:
                info[4] = start

            if oldEnd < end:
                info[5] = end

            exonLst.append((start, end))
        else:
            exonLst = [(start, end)]
            ujcDct[xscript] = [exonLst,
                               xscript,
                               geneID,
                               seqname,
                               start,
                               end,
                               strand,
                               source]

    for x, info in ujcDct.items():
        startValues, endValues = zip(*sorted(info[0]))
        jxns = list(zip(endValues[:-1], startValues[1:]))

        if jxns == []:
            jxnStr = "{}_{}_monoexon_{}_{}".format(
                info[3], info[6], info[4], info[5])

            numJxn = 0
        else:
            jxnLst = []
            for jxn in jxns:
                jxnLst.append("{}_{}".format(jxn[0], jxn[1]))
            jxnStr = "_".join(jxnLst)

            jxnStr = "{}_{}_".format(info[3], info[6]) + jxnStr

            numJxn = len(jxns)

        info.append(jxnStr)
        info.append(numJxn)

    return ujcDct


def createUJCIndex(ujcDct, keepGene=False, trackSrc=False):
    """
    Takes extracted junction information and creates a dataframe that is 
    UJC focused (all transcripts under one UJC grouped into the transcript_id column).

    Parameters
    ----------
    ujcDct : DICTIONARY {Transcript_id: [info]}
            A dictionary of transcripts keyed to their info.
    trPrefix: STRING
            A prefix for the ujc_id.

    Returns
    -------
    allUJC : DATAFRAME
            Dataframe with information on the UJCs, with their ids, transcripts, etc.
    """

    monoExonDct = dict()
    multiExonDct = dict()

    for xscript, info in ujcDct.items():
        jxnStr = info[8]

        if 'monoexon' not in jxnStr:
            multiExonDct.update({xscript: info})
        else:
            monoExonDct.update({xscript: info})

    if len(monoExonDct) > 0:
        monoXscriptDf = pd.DataFrame(
            monoExonDct,
            index=pd.Index(["exons",
                            "transcriptID",
                            "geneID",
                            "chr",
                            "start",
                            "end",
                            "strand",
                            "source",
                            "jxnString",
                            "numJxn"])
        ).T.sort_values(by=["jxnString", "start", "end"])

        # THIS IS WHERE OVERLAPPING MONOEXONS IN THE SAME GENE ARE COLLAPSED
        monoXscriptDf['tmpStart'] = monoXscriptDf['start']
        monoXscriptDf['tmpEnd'] = monoXscriptDf['end']

        appendedRowLst = []
        for row in monoXscriptDf.to_dict('records'):
            if appendedRowLst:
                lastRow = appendedRowLst[-1]

                if lastRow['chr'] == row['chr'] and lastRow['strand'] == row['strand'] and lastRow['geneID'] == row['geneID']:
                    if lastRow['tmpEnd'] > row['tmpStart']:

                        row['tmpStart'] = lastRow['tmpStart']

                        if (lastRow['tmpEnd'] < row['tmpEnd']):
                            for loopRow in appendedRowLst:
                                if loopRow['chr'] == row['chr'] and loopRow['strand'] == row['strand'] and lastRow['geneID'] == row['geneID']:
                                    loopRow['tmpEnd'] = row['tmpEnd']
                        else:
                            row['tmpEnd'] = lastRow['tmpEnd']

                        appendedRowLst.append(row)
                    else:
                        appendedRowLst.append(row)
                else:
                    appendedRowLst.append(row)
            else:
                appendedRowLst.append(row)

        for row in appendedRowLst:
            jString = ("{}_{}_monoexon_{}_{}".format(row['chr'],
                                                     row['strand'],
                                                     row['tmpStart'],
                                                     row['tmpEnd']))

            row['jxnString'] = jString

        newMonoDf = pd.DataFrame(appendedRowLst)

        monoUJC = newMonoDf.sort_values(by=['start', 'end'])
        monoUJC.drop(columns=['tmpStart', 'tmpEnd'])

        monoUJC['pair'] = list(zip(monoUJC['geneID'], monoUJC['transcriptID']))

        monoUJC = monoUJC.groupby(["jxnString"]).agg({
            "geneID": set,
            "chr": "first",
            "start": "min",
            "end": "max",
            "strand": "first",
            "numJxn": "max",
            "transcriptID": set,
            "pair": lambda x: list(set(x))}).reset_index()

    else:
        monoExonDct = None

    if len(multiExonDct) > 0:

        multiXscriptDf = pd.DataFrame(
            multiExonDct,
            index=pd.Index(["exons",
                            "transcriptID",
                            "geneID",
                            "chr",
                            "start",
                            "end",
                            "strand",
                            "source",
                            "jxnString",
                            "numJxn"])
        ).T.sort_values(by=["jxnString", "start", "end"])

        multiXscriptDf['pair'] = list(
            zip(multiXscriptDf['geneID'], multiXscriptDf['transcriptID']))

        # THIS IS WHERE TRANSCRIPTS WITH 1+ EXONS HAVE THEIR STARTS/ENDS COLLAPSED
        multiUJC = multiXscriptDf.groupby(["jxnString"]).agg({
            "geneID": set,
            "chr": "first",
            "start": "min",
            "end": "max",
            "strand": "first",
            "numJxn": "max",
            "transcriptID": set,
            "pair": lambda x: list(set(x))}).reset_index()

    else:
        multiExonDct = None

    if monoExonDct and multiExonDct:
        ujcDscrptnDf = pd.concat([monoUJC, multiUJC], ignore_index=True)
    elif monoExonDct:
        ujcDscrptnDf = monoUJC.copy()
        del (monoUJC)
    else:
        ujcDscrptnDf = multiUJC.copy()
        del (multiUJC)

    ujcDscrptnDf['jxnHash'] = ujcDscrptnDf['jxnString'].apply(
        lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest())

    print("Number of UJCs: {}".format(len(ujcDscrptnDf['jxnHash'])))

    if not ujcDscrptnDf['jxnHash'].is_unique:
        print("Wow! A rare jxnHash collision: two jxnStrings have resulted in the exact same hash for these genes and transcripts: ")
        print("geneID", "transcriptID")

        duplicateDf = ujcDscrptnDf[ujcDscrptnDf.duplicated(
            subset='jxnHash', keep=False) | ujcDscrptnDf.duplicated(subset='jxnHash', keep='first')]
        for row in duplicateDf.to_dict('records'):
            print(row['geneID'], row['transcriptID'])

    ujcDscrptnDf['flagMultiGene'] = ujcDscrptnDf['geneID'].apply(
        lambda x: 1 if len(x) > 1 else 0)
    ujcDscrptnDf['flagMultiXscript'] = ujcDscrptnDf['pair'].apply(lambda pair: 1 if len(
        set([tup[0] for tup in pair])) < len([tup[0] for tup in pair]) else 0)

    ujcDscrptnDf = ujcDscrptnDf.sort_values(
        by=['chr', 'strand', 'start'], ascending=True)

    if keepGene:
        ujcOutDf = ujcDscrptnDf[['jxnHash', 'geneID', 'flagMultiXscript', 'flagMultiGene',
                                 'numJxn', 'chr', 'strand', 'start', 'end', 'jxnString']].copy()

        ujcOutDf['geneID'] = ujcOutDf['geneID'].apply(lambda x: '|'.join(x))
    else:
        ujcOutDf = ujcDscrptnDf[['jxnHash', 'flagMultiXscript', 'flagMultiGene',
                                 'numJxn', 'chr', 'strand', 'start', 'end', 'jxnString']]

    ujcOutDf = ujcOutDf.rename(
        columns={'start': 'donorStart', 'end': 'acceptorEnd'})

    xscriptLinkDf = ujcDscrptnDf.copy(deep=True).reset_index()
    xscriptLinkDf = xscriptLinkDf[['pair', 'jxnHash', 'jxnString']]
    xscriptLinkDf = xscriptLinkDf.explode('pair')

    # Left for testing later (much faster)
    # xscriptLinkDf[['geneID', 'transcriptID']] = pd.DataFrame(
    #     xscriptLinkDf['pair'].tolist(), index=xscriptLinkDf.index)

    xscriptLinkDf[['geneID', 'transcriptID']
                  ] = xscriptLinkDf['pair'].apply(pd.Series)

    xscriptLinkDf = xscriptLinkDf.drop_duplicates(
    )[['transcriptID', 'geneID', 'jxnHash', 'jxnString']]

    if trackSrc:
        xscript2SrcDf = pd.DataFrame([
            (xscript, info[7]) for xscript, info in ujcDct.items()
        ], columns=['transcriptID', 'source'])

        xscriptLinkDf = pd.merge(xscriptLinkDf, xscript2SrcDf,
                                 on='transcriptID', how='outer', indicator='merge_check')

        if not (xscriptLinkDf['merge_check'] == "both").all():
            raise Exception(
                "An error occurred when adding sample column to xscript_link file.")
        else:
            xscriptLinkDf = xscriptLinkDf.drop('merge_check', axis=1)

        xscriptLinkDf = xscriptLinkDf.drop_duplicates(
        )[['source', 'transcriptID', 'geneID', 'jxnHash', 'jxnString']]

    return ujcDscrptnDf, ujcOutDf, xscriptLinkDf


def createExonOutput(ujcDf, ujcDct, keepGene=False):
    """
    Creates the dataframe with exon information to be output as a GTF file
    using the UJCs as transcripts.

    Parameters
    ----------
    ujcDf : DATAFRAME
            Dataframe with information on the UJCs, with their ids, transcripts, etc.

    ujcDct : DICTIONARY {Transcript_id: [info]}
            A dictionary of transcripts keyed to their info.

    Returns
    -------
    outExonDf : DATAFRAME
            A dataframe in the proper format to be written as a GTF file.

    """

    workingDf = ujcDf.explode(
        'pair')[['pair', 'chr', 'strand', 'jxnHash', 'start', 'end', 'numJxn']]
    workingDf[['geneID', 'transcriptID']] = pd.DataFrame(
        workingDf['pair'].to_list(), index=workingDf.index)
    workingDf.drop('pair', axis=1, inplace=True)

    if keepGene:
        workingDf = workingDf.groupby(['jxnHash']).agg({
            "chr": "first",
            "start": "min",
            "end": "max",
            "strand": "first",
            "transcriptID": set,
            "geneID": set,
            "numJxn": "max"}).reset_index()

        if workingDf['geneID'].apply(lambda x: len(x) > 1).any():
            raise Exception(
                "Error: The keepGene parameter is on but the GTF contains more than one gene.")
        else:
            workingDf['geneID'] = workingDf['geneID'].apply(
                lambda x: next(iter(x)))

    else:
        workingDf = workingDf.groupby(['jxnHash']).agg({
            "chr": "first",
            "start": "min",
            "end": "max",
            "strand": "first",
            "transcriptID": set,
            "numJxn": "max"}).reset_index()

    seqnameLst = []
    startLst = []
    endLst = []
    strandLst = []
    hashLst = []
    geneIDLst = []

    # tested -> ujcDf contains accurate start and end
    for row in workingDf.to_dict('records'):

        seqname = row['chr']
        strand = row['strand']
        jxnHash = row['jxnHash']

        if keepGene:
            geneID = row['geneID']
        else:
            geneID = row['jxnHash']

        firstStart = row['start']
        lastEnd = row['end']

        # tested and proved all xscripts under same UJC have same junctions and internal exons

        xscript = next(iter(row['transcriptID']))

        exons = ujcDct[xscript][0]

        flagMono = row['numJxn'] < 1

        if flagMono:
            seqnameLst.append(seqname)
            startLst.append(firstStart)
            endLst.append(lastEnd)
            hashLst.append(jxnHash)
            strandLst.append(strand)
            geneIDLst.append(geneID)
        else:
            seqnameLst.append(seqname)
            startLst.append(firstStart)
            endLst.append(exons[0][1])
            hashLst.append(jxnHash)
            strandLst.append(strand)
            geneIDLst.append(geneID)

            for exon in exons[1:-1]:
                seqnameLst.append(seqname)
                startLst.append(exon[0])
                endLst.append(exon[1])
                hashLst.append(jxnHash)
                strandLst.append(strand)
                geneIDLst.append(geneID)

            seqnameLst.append(seqname)
            startLst.append(exons[-1][0])
            endLst.append(lastEnd)
            hashLst.append(jxnHash)
            strandLst.append(strand)
            geneIDLst.append(geneID)

    outExonDf = pd.DataFrame(
        {
            'seqname': seqnameLst,
            'start': startLst,
            'end': endLst,
            'strand': strandLst,
            'transcript_id': hashLst,
            'gene_id': geneIDLst
        })

    outExonDf = outExonDf.sort_values(by=['seqname', 'transcript_id', 'start'])

    # numColumns = ['start', 'end']
    # outExonDf[numColumns] = outExonDf[numColumns].astype(int)
    # result = outExonDf[outExonDf['end'] < outExonDf['start']]

    return outExonDf


def main():
    """
    Run the program.

    Returns
    -------
    None.

    """

    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dmel_fb650/dmel-all-r6.50.gtf"
    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_id_ujc_update/subset_dm650.gtf"
    # outdir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_id_ujc_update"
    # includeGTF = True
    # includeCnt = True
    # keepGene = True
    # trackSrc = True

    inGTF = args.inGTF
    outdir = args.outdir
    includeCnt = args.includeCnt
    includeGTF = args.includeGTF
    keepGene = args.keepGene
    trackSrc = args.trackSrc

    print("Loading...")
    alphatic = time.perf_counter()

    exonData = trand.io.read_exon_data_from_file(
        infile=inGTF, keepSrc=trackSrc)
    prefix = os.path.basename(inGTF).split('.')[0]

    exonData = exonData[exonData['gene_id'] == "FBgn0264270"]

    toc = time.perf_counter()
    print(
        f"GTF Read complete! Took {toc-alphatic:0.4f} seconds. Extracting junctions...")
    tic = time.perf_counter()

    ujcDct = extractJunction(exonData=exonData, trackSrc=trackSrc)

    toc = time.perf_counter()
    print(
        f"Complete! Operation took {toc-tic:0.4f} seconds. Creating UJC DataFrame...")
    tic = time.perf_counter()

    ujcDf, dscDf, linkDf = createUJCIndex(
        ujcDct=ujcDct, keepGene=keepGene, trackSrc=trackSrc)

    toc = time.perf_counter()
    print(f"Complete! Operation took {toc-tic:0.4f} seconds. Writing files...")
    tic = time.perf_counter()

    dscOutPath = outdir + "/" + prefix + "_ujc_dscrptn.csv"
    linkOutPath = outdir + "/" + prefix + "_ujc_xscript_link.csv"

    try:
        dscDf.to_csv(dscOutPath, index=False)
        linkDf.to_csv(linkOutPath, index=False)
    except OSError:
        raise OSError("Output directory must already exist.")

    if includeGTF:

        print("Writing GTF...")

        gtfDf = createExonOutput(ujcDf=ujcDf, ujcDct=ujcDct, keepGene=keepGene)
        gtfOutPath = outdir + "/" + prefix + "_ujc.gtf"

        if os.path.isfile(gtfOutPath):
            os.remove(gtfOutPath)

        trand.io.write_gtf(data=gtfDf, out_fhs={
                           "gtf": gtfOutPath}, fh_name="gtf")

    if includeCnt:

        print("Counting transcripts per UJC...")

        if trackSrc:
            countDf = linkDf.groupby(['jxnHash', 'source']).count()[
                'transcriptID'].reset_index()

            countDf = countDf.rename(
                columns={'transcriptID': 'numTranscripts'})
            countDf = countDf[['source', 'jxnHash', 'numTranscripts']]

        else:
            countDf = linkDf.groupby(['jxnHash']).count()[
                'transcriptID'].reset_index()

            countDf.columns = ['jxnHash', 'numTranscripts']

        countOutPath = outdir + "/" + prefix + "_ujc_count.csv"

        try:
            countDf.to_csv(countOutPath, index=False)
        except OSError:
            raise OSError("Output directory must already exist.")

    omegatoc = time.perf_counter()
    print(f"Complete! Operation took {omegatoc-alphatic:0.4f} total seconds.")


if __name__ == '__main__':
    global args
    args = getOptions()
    main()

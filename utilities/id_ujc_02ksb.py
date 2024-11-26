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


# TODO:




import trand.io
import sys
import copy
import hashlib
import csv
import os
import pandas as pd
import time
import argparse
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
                                     "1. An info file that is unique on jxnHash and provides the following "
                                     "information: jxnHash, flagMultiXscript, flagMultiGene, chr, strand, numJxn, donorStart, acceptorEnd, jxnString "
                                     "and jxnHash. jxnStrings are in the format of chr_strand_start:end_start:end_start:end "
                                     "with monoexons as chr_strand_monoexon_start_end (with the start/end being for the exon)."
                                     "jxnHash is the unique identifier for the group of transcripts. 2. An xscript link"
                                     "file unique on transcriptID that more directly links a read/transcript to its "
                                     "jxnHash. 3. A GTF file with representative transcript models "
                                     "for each group. The jxnHash will be the \'transcript_id\' for the group. "
                                     "Input: a GTF file (--gtf), and an output directory (--outdir). Output files will begin with "
                                     "the name of the input GTF file. Allows the option to skip the output of the GTF file "
                                     "with representative transcript models. (--skip-gtf). Allows the option to output "
                                     "another key file with the number of transcripts per jxnHash counted (--count-ujc). "
                                     "Add -kg (--keep-geneID) to include geneID in output. There can jxnHash with more than one gene. "
                                     "This option is most useful when running id_ujc on a GTF split by gene. "
                                     "Add -ts (--track-source) to include the source column of the GTF in the output. "
                                     "This option is most useful if the source column contains information about which sample "
                                     "each transcript orignated from.")

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

    exonDfr = checkStrandAndChromosome(exonData=exonData)

    print("Number of transcripts: ", end="")
    print(len(exonDfr['transcript_id'].unique()))

    print("Number of genes: ", end="")
    print(len(exonDfr['gene_id'].unique()))

    # First, instead of grouping, then sorting
    # Sort by transcript -> sort by start. the whole dataframe
    sortedDfr = exonDfr.sort_values(
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

    for row in sortedDfr.to_dict('records'):
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


def createUJCCSV(ujcDct, keepGene=False, trackSrc=False):
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
        monoXscriptDfr = pd.DataFrame(
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
        monoXscriptDfr['tmpStart'] = monoXscriptDfr['start']
        monoXscriptDfr['tmpEnd'] = monoXscriptDfr['end']

        appendedRowLst = []
        for row in monoXscriptDfr.to_dict('records'):
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

        newMonoDfr = pd.DataFrame(appendedRowLst)

        monoUJC = newMonoDfr.sort_values(by=['start', 'end'])
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

        multiXscriptDfr = pd.DataFrame(
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

        multiXscriptDfr['pair'] = list(
            zip(multiXscriptDfr['geneID'], multiXscriptDfr['transcriptID']))

        # THIS IS WHERE TRANSCRIPTS WITH 1+ EXONS HAVE THEIR STARTS/ENDS COLLAPSED
        multiUJC = multiXscriptDfr.groupby(["jxnString"]).agg({
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
        ujcInfoDfr = pd.concat([monoUJC, multiUJC], ignore_index=True)
    elif monoExonDct:
        ujcInfoDfr = monoUJC.copy()
        del (monoUJC)
    else:
        ujcInfoDfr = multiUJC.copy()
        del (multiUJC)

    ujcInfoDfr['jxnHash'] = ujcInfoDfr['jxnString'].apply(
        lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest())

    print("Number of UJCs: {}".format(len(ujcInfoDfr['jxnHash'])))

    if not ujcInfoDfr['jxnHash'].is_unique:
        print("Wow! A rare jxnHash collision: two jxnStrings have resulted in the exact same hash for these genes and transcripts: ")
        print("geneID", "transcriptID")

        duplicateDfr = ujcInfoDfr[ujcInfoDfr.duplicated(
            subset='jxnHash', keep=False) | ujcInfoDfr.duplicated(subset='jxnHash', keep='first')]
        for row in duplicateDfr.to_dict('records'):
            print(row['geneID'], row['transcriptID'])

    ujcInfoDfr['flagMultiGene'] = ujcInfoDfr['geneID'].apply(
        lambda x: 1 if len(x) > 1 else 0)
    ujcInfoDfr['flagMultiXscript'] = ujcInfoDfr['pair'].apply(lambda pair: 1 if len(
        set([tup[0] for tup in pair])) < len([tup[0] for tup in pair]) else 0)

    ujcInfoDfr = ujcInfoDfr.sort_values(
        by=['chr', 'strand', 'start'], ascending=True)

    if keepGene:
        ujcOutDfr = ujcInfoDfr[['jxnHash', 'geneID', 'flagMultiXscript', 'flagMultiGene',
                                'numJxn', 'chr', 'strand', 'start', 'end', 'jxnString']].copy()

        ujcOutDfr['geneID'] = ujcOutDfr['geneID'].apply(lambda x: '|'.join(x))
    else:
        ujcOutDfr = ujcInfoDfr[['jxnHash', 'flagMultiXscript', 'flagMultiGene',
                                'numJxn', 'chr', 'strand', 'start', 'end', 'jxnString']]

    ujcOutDfr = ujcOutDfr.rename(
        columns={'start': 'donorStart', 'end': 'acceptorEnd'})

    xscriptLinkDfr = ujcInfoDfr.copy(deep=True).reset_index()
    xscriptLinkDfr = xscriptLinkDfr[['pair', 'jxnHash', 'jxnString']]
    xscriptLinkDfr = xscriptLinkDfr.explode('pair')

    # Left for testing later (much faster)
    # xscriptLinkDfr[['geneID', 'transcriptID']] = pd.DataFrame(
    #     xscriptLinkDfr['pair'].tolist(), index=xscriptLinkDfr.index)

    xscriptLinkDfr[['geneID', 'transcriptID']
                   ] = xscriptLinkDfr['pair'].apply(pd.Series)

    xscriptLinkDfr = xscriptLinkDfr.drop_duplicates(
    )[['transcriptID', 'geneID', 'jxnHash']]

    if trackSrc:
        xscript2SrcDfr = pd.DataFrame([
            (xscript, info[7]) for xscript, info in ujcDct.items()
        ], columns=['transcriptID', 'source'])

        xscriptLinkDfr = pd.merge(xscriptLinkDfr, xscript2SrcDfr,
                                  on='transcriptID', how='outer', indicator='merge_check')

        if not (xscriptLinkDfr['merge_check'] == "both").all():
            raise Exception(
                "An error occurred when adding sample column to xscript_link file.")
        else:
            xscriptLinkDfr = xscriptLinkDfr.drop('merge_check', axis=1)

        xscriptLinkDfr = xscriptLinkDfr.drop_duplicates(
        )[['source', 'transcriptID', 'geneID', 'jxnHash']]

    return ujcInfoDfr, ujcOutDfr, xscriptLinkDfr


def createExonOutput(infoDfr, ujcDct, keepGene=False):
    """
    Creates the dataframe with exon information to be output as a GTF file
    using the UJCs as transcripts.

    Parameters
    ----------
    ujcDfr : DATAFRAME
            Dataframe with information on the UJCs, with their ids, transcripts, etc.

    ujcDct : DICTIONARY {Transcript_id: [info]}
            A dictionary of transcripts keyed to their info.

    Returns
    -------
    outExonDfr : DATAFRAME
            A dataframe in the proper format to be written as a GTF file.

    """

    workingDfr = infoDfr.explode(
        'pair')[['pair', 'chr', 'strand', 'jxnHash', 'start', 'end', 'numJxn']]
    workingDfr[['geneID', 'transcriptID']] = pd.DataFrame(
        workingDfr['pair'].to_list(), index=workingDfr.index)
    workingDfr.drop('pair', axis=1, inplace=True)

    if keepGene:
        workingDfr = workingDfr.groupby(['jxnHash']).agg({
            "chr": "first",
            "start": "min",
            "end": "max",
            "strand": "first",
            "transcriptID": set,
            "geneID": set,
            "numJxn": "max"}).reset_index()

        if workingDfr['geneID'].apply(lambda x: len(x) > 1).any():
            raise Exception(
                "Error: The keepGene parameter is on but the GTF contains more than one gene.")
        else:
            workingDfr['geneID'] = workingDfr['geneID'].apply(
                lambda x: next(iter(x)))

    else:
        workingDfr = workingDfr.groupby(['jxnHash']).agg({
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

    # tested -> ujcDfr contains accurate start and end
    for row in workingDfr.to_dict('records'):

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

    outExonDfr = pd.DataFrame(
        {
            'seqname': seqnameLst,
            'start': startLst,
            'end': endLst,
            'strand': strandLst,
            'transcript_id': hashLst,
            'gene_id': geneIDLst
        })

    outExonDfr = outExonDfr.sort_values(
        by=['seqname', 'transcript_id', 'start'])

    # numColumns = ['start', 'end']
    # outExonDfr[numColumns] = outExonDfr[numColumns].astype(int)
    # result = outExonDfr[outExonDfr['end'] < outExonDfr['start']]

    return outExonDfr


def main():
    """
    Run the program.

    Returns
    -------
    None.

    """

    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/dmel_fb650/dmel-all-r6.50.gtf"
    # # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_id_ujc_update/subset_dm650.gtf"
    # outdir = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_id_ujc_update"
    # includeGTF = True
    # includeCnt = True
    # keepGene = True
    # trackSrc = True

    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/transcript_ortholog/ujc_species_gtf/roz_dyak_data_XLOC_018390/dyak_data_XLOC_018390_job__run_26080.gtf"

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

    toc = time.perf_counter()
    print(
        f"GTF Read complete! Took {toc-alphatic:0.4f} seconds. Extracting junctions...")
    tic = time.perf_counter()

    ujcDct = extractJunction(exonData=exonData, trackSrc=trackSrc)

    toc = time.perf_counter()
    print(
        f"Complete! Operation took {toc-tic:0.4f} seconds. Creating UJC DataFrame...")
    tic = time.perf_counter()

    infoDfr, outInfoDfr, linkDfr = createUJCCSV(
        ujcDct=ujcDct, keepGene=keepGene, trackSrc=trackSrc)

    toc = time.perf_counter()
    print(f"Complete! Operation took {toc-tic:0.4f} seconds. Writing files...")
    tic = time.perf_counter()

    infoOutPath = outdir + "/" + prefix + "_ujc_info.csv"
    linkOutPath = outdir + "/" + prefix + "_ujc_xscript_link.csv"

    try:
        outInfoDfr.to_csv(infoOutPath, index=False)
        linkDfr.to_csv(linkOutPath, index=False)
    except OSError:
        raise OSError("Output directory must already exist.")

    if includeGTF:

        print("Writing GTF...")

        gtfDfr = createExonOutput(
            infoDfr=infoDfr, ujcDct=ujcDct, keepGene=keepGene)
        gtfOutPath = outdir + "/" + prefix + "_ujc.gtf"

        if os.path.isfile(gtfOutPath):
            os.remove(gtfOutPath)

        trand.io.write_gtf(data=gtfDfr, out_fhs={
                           "gtf": gtfOutPath}, fh_name="gtf")

    if includeCnt:

        print("Counting transcripts per UJC...")

        if trackSrc:
            countDfr = linkDfr.groupby(['jxnHash', 'source']).count()[
                'transcriptID'].reset_index()

            countDfr = countDfr.rename(
                columns={'transcriptID': 'numTranscripts'})
            countDfr = countDfr[['source', 'jxnHash', 'numTranscripts']]

        else:
            countDfr = linkDfr.groupby(['jxnHash']).count()[
                'transcriptID'].reset_index()

            countDfr.columns = ['jxnHash', 'numTranscripts']

        countOutPath = outdir + "/" + prefix + "_ujc_count.csv"

        try:
            countDfr.to_csv(countOutPath, index=False)
        except OSError:
            raise OSError("Output directory must already exist.")

    omegatoc = time.perf_counter()
    print(f"Complete! Operation took {omegatoc-alphatic:0.4f} total seconds.")


if __name__ == '__main__':
    global args
    args = getOptions()
    main()

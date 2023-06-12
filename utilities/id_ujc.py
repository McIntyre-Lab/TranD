#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on Fri May 29 13:34:55 2023

@author: k.bankole
"""

"""

Identify the unique junction chains (UJC) of a GTF file and combine
transcripts into UJCs.

Created from a previous utility in TranD named consolidation

Version 9: Change the ujc_id to be a representative transcript.

"""

import argparse
import time
import pandas as pd
import os
import pickle
import trand.io
from dataclasses import dataclass

@dataclass()
class J_CHAIN:
        """
        Object representation of a junction string.
        Makes it easier to extract junctions and provides ability output them
        as a string, just like before (chainToStr).        
        
        seqname: seqname
        junctionLst: a list of tuples (lastExonEnd, nextExonStart) for each junction
        strand: strand
        """
        
        seqname: str
        junctionLst: list
        strand: str
        
        # junctionLst format, for posterity:
                # (lastExonEnd, nextExonStart)
                # just a list of pairs of values like this
        
        
        def chainToStr(self):
                chainStr = []
                
                for junctionPair in self.junctionLst:
                        chainStr.append(self.seqname + ":" + str(junctionPair[0]) + ":" + str(junctionPair[1])
                                        + ":" + self.strand)
                        
                return "|".join(chainStr)

        

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
                                         "GTF file and combines transcripts based on these UJCs. Outputs a summary file "
                                         "containing info on these combined transcripts and their \"ujc_id.\" "
                                         "Please note: The ujc_id is a representative transcript for all the "
                                         "transcripts with that same junction chain. The utility sorts the group of "
                                         "transcripts alphabetically and selects the first one as the representative."
                                         "Also, outputs a GTF file with representative transcript models for each UJC. Input a GTF "
                                         "file (--gtf), an output directory (--outdir) and a prefix for the output files (--prefix). "
                                         "Allows the option to skip the output of the GTF file with representative transcript models."
                                         "(--skip-gtf).")
        
        ## INPUT
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
                dest="outGTF",
                action="store_false",
                help="Use this argument to remove the output of a GTF with "
                "representative transcript models for each UJC."
                        "Defaults to outputting the GTF."
        )
        
        ## OUTPUT
        parser.add_argument(
                "-o",
                "--outdir",
                dest="outdir",
                required=True,
                help="Location of output directory, must already exist."
        )
        
        parser.add_argument(
                "-p",
                "--prefix",
                dest="prefix",
                required=True,
                help="Prefix for the output file(s). Example: prefix_UJC_ID.csv"
        )
        
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
                badStrand = list(strandCheck[strandCheck>1].index)
                for gene in badStrand:
                        print("!!! WARNING: gene {} contains transcripts/exons on both strands - "
                              "removing.".format(gene))
                        
                exonData = exonData[~exonData["gene_id"].isin(badStrand)]
        
        
        chrCheck = geneGrps["seqname"].nunique()
        if (chrCheck > 1).any():
            badChr = list(chrCheck[chrCheck>1].index)
            for gene in badChr:
                    print("!!! WARNING: gene {} contains transcripts/exons on difference chromosomes - "
                          "removing.".format(gene))
                    
            exonData = exonData[~exonData["gene_id"].isin(badChr)]
            
        return exonData


def extractJunction(exonData):
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
        
        print ("Number of transcripts: ", end="")
        print (len(exonDf['transcript_id'].unique()))
        
        # First, instead of grouping, then sorting
        # Sort by transcript -> sort by start. the whole dataframe 
        sortedDf = exonDf.sort_values(by=['transcript_id', 'start']).reset_index(drop=True)
                        
        ujcDct = {}
        # Value Legend:
                # 0 = exons
                # 1 = xscript
                # 2 = gene
                # 3 = seqname
                # 4 = start
                # 5 = end
                # 6 = strand
        
        for row in sortedDf.to_dict('records'):                
                xscript = row['transcript_id']
                
                seqname = row['seqname']
                strand = row['strand']
                geneID = row['gene_id']
                
                start = row['start']
                end = row['end']
                
                
                if xscript in ujcDct.keys():
                        info = ujcDct[xscript]
                        exonLst = info[0]
                        oldStart = info[4]
                        oldEnd = info[5]
                        
                        if oldStart > start:
                                info[4] = start
                        
                        if oldEnd < end:
                                info[5] = end
                        
                        exonLst.append((start,end))
                else:
                        exonLst = [(start,end)]
                        ujcDct[xscript] = [exonLst,
                                            xscript,
                                            geneID,
                                            seqname,
                                            start,
                                            end,
                                            strand]        
        
        
        
        for x, info in ujcDct.items():
                startValues, endValues = zip(*sorted(info[0]))
                junctions = list(zip(endValues[:-1], startValues[1:]))
                
                if junctions == []:
                        junctionChain = None
                else:    
                        junctionChain = J_CHAIN(seqname=info[3], junctionLst=junctions, strand=info[6])
                info.append(junctionChain)
        
        return ujcDct

def createUJCDf(ujcDct):
        """
        Takes extracted junction information and creates a dataframe that is 
        UJC focused (all transcripts under one UJC grouped into the transcript_id column).

        Parameters
        ----------
        ujcDct : DICTIONARY {Transcript_id: [info]}
                A dictionary of transcripts keyed to their info.

        Returns
        -------
        allUJC : DATAFRAME
                Dataframe with information on the UJCs, with their ids, transcripts, etc.
        """
        
        
        monoExonDct = dict()
        multiExonDct = dict()        
        
        for xscript, info in ujcDct.items():
                junctionChain = info[7]
                
                if junctionChain:
                        newInfo = [junctionChain.chainToStr(), info[1], info[2], info[3], info[4], info[5], info[6]]
                        multiExonDct.update({xscript: newInfo})
                else:
                        newInfo = [info[0], info[1], info[2], info[3], info[4], info[5], info[6]]
                        monoExonDct.update({xscript: newInfo})
                
        
        if len (monoExonDct) > 0:
                monoXscriptDf = pd.DataFrame(monoExonDct,
                                                index = pd.Index(["junction_string",
                                                                  "transcript_id",
                                                                  "gene_id",
                                                                  "seqname",
                                                                  "start", "end", "strand"])
                                                ).T.sort_values(by=["gene_id", "start", "end"])
                                                
                monoXscriptDf['tmpStart'] = monoXscriptDf['start']
                monoXscriptDf['tmpEnd'] = monoXscriptDf['end']
                
                appendedRowLst = []
                for row in monoXscriptDf.to_dict('records'):
                        if appendedRowLst:
                                lastRow = appendedRowLst[-1]
                                        
                                if lastRow['gene_id'] == row['gene_id']:
                                        if lastRow['tmpEnd'] > row['tmpStart']:
                                                
                                                row['tmpStart'] = lastRow['tmpStart']
                                                
                                                if (lastRow['tmpEnd'] < row['tmpEnd']):
                                                        for loopRow in appendedRowLst:
                                                                if loopRow['gene_id'] == row['gene_id']:
                                                                        loopRow['tmpEnd'] = row['tmpEnd']
                                                else:
                                                        row['tmpEnd'] = lastRow['tmpEnd']
                                                        
                                                appendedRowLst.append(row)
                                else:
                                        appendedRowLst.append(row)
                        else:
                                appendedRowLst.append(row)
                                
                for row in appendedRowLst:
                        jString = ("monoexon_"
                           + str(row['tmpStart']) + "_"
                           + str(row['tmpEnd']))
                        
                        row['junction_string'] = jString
                
                newMonoDf = pd.DataFrame(appendedRowLst)
                
                monoUJC = newMonoDf.sort_values(by=['start', 'end'])
                monoUJC.drop(columns=['tmpStart','tmpEnd'])
                monoUJC = monoUJC.groupby(["gene_id", "junction_string"]).agg({
                        "seqname":"first",
                        "start":"min",
                        "end":"max",
                        "strand":"first",
                        "transcript_id": lambda x: '|'.join(x)}).reset_index()                
        else:
                monoExonDct = None
        
        if len(multiExonDct) > 0:
                multiUJC = pd.DataFrame(multiExonDct, index=pd.Index(["junction_string", 
                                                                    "transcript_id", 
                                                                    "gene_id", 
                                                                    "seqname", 
                                                                    "start", "end", "strand"])
                                        ).T.groupby(["gene_id","junction_string"]).agg({
                        "seqname": "first",
                        "start": "min",
                        "end": "max",
                        "strand": "first",
                        "transcript_id": lambda x: "|".join(x)}).reset_index()
                                
        else:
                multiExonDct = None
                
        if monoExonDct and multiExonDct:
                allUJC = pd.concat([monoUJC, multiUJC], ignore_index=True)
        elif monoExonDct:
                allUJC = monoUJC.copy()
                del(monoUJC)
        else:
                allUJC = multiUJC.copy()
                del(multiUJC)
                
                
        sort_order = {"gene_id": "asc", "start": "asc", "transcript_id": "asc"}
        allUJC = allUJC.sort_values(by=list(sort_order.keys()), ascending=[True if val=="asc" else False for val in sort_order.values()])
                
        allUJC["split"] = allUJC["transcript_id"].str.split('|').apply(sorted)
        
        allUJC["ujc_id"] = allUJC["split"].str[0]
        
        allUJC.drop(columns="split")
                
        return allUJC

def createExonOutput(ujcDf, ujcDct):
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
        exonDf : DATAFRAME
                A dataframe in the proper format to be written as a GTF file.

        """
        seqnameLst = []
        startLst = []
        endLst = []
        strandLst = []
        ujcIDLst = []
        geneIDLst = []
        
        # tested -> ujcDf contains accurate start and end        
        for row in ujcDf.to_dict('records'):
                seqname = row['seqname']
                strand = row['strand']
                ujcID = row['ujc_id']
                geneID = row['gene_id']
                
                firstStart = row['start']
                lastEnd = row['end']
                
                # tested and proved all xscripts under same UJC have same junctions and internal exons
                xscript = row['transcript_id'].split('|')[0]
                
                exons = ujcDct[xscript][0]
                flagMono =  ujcDct[xscript][7] == None
                

                
                if flagMono:
                        seqnameLst.append(seqname)
                        startLst.append(firstStart)
                        endLst.append(lastEnd)
                        ujcIDLst.append(ujcID)
                        strandLst.append(strand)
                        geneIDLst.append(geneID)
                else:
                        startValues, endValues = zip(*sorted(exons))
                        
                        seqnameLst.append(seqname)
                        startLst.append(firstStart)
                        endLst.append(endValues[0])
                        ujcIDLst.append(ujcID)
                        strandLst.append(strand)
                        geneIDLst.append(geneID)
                        
                        for idx in range(1, len(exons)-2):
                              seqnameLst.append(seqname)
                              startLst.append(startValues[idx])
                              endLst.append(endValues[idx])
                              ujcIDLst.append(ujcID)
                              strandLst.append(strand)
                              geneIDLst.append(geneID)      
                        
                        seqnameLst.append(seqname)
                        startLst.append(startValues[len(exons)-1])
                        endLst.append(lastEnd)
                        ujcIDLst.append(ujcID)
                        strandLst.append(strand)
                        geneIDLst.append(geneID)
                
        exonDf = pd.DataFrame(
                {
                        'seqname':seqnameLst,
                        'start':startLst,
                        'end':endLst,
                        'strand':strandLst,
                        'transcript_id':ujcIDLst,
                        'gene_id':geneIDLst
                })
                
        return exonDf

def createOutput(ujcDf, ujcDct):
        """

        Parameters
        ----------
        ujcDf : DATAFRAME
                Dataframe with information on the UJCs, with their ids, transcripts, etc.
                
        ujcDct : DICTIONARY {Transcript_id: [info]}
                A dictionary of transcripts keyed to their info.

        Returns
        -------
        outDf : DATAFRAME
                A dataframe that is the list of transcripts and which UJCs those transcripts
                belong to.
        
        """
        xscriptLst = []
        ujcIDLst = []
        geneLst = []
        junctionLst = []
        
        for row in ujcDf.to_dict('records'):
                xscripts = row['transcript_id'].split('|')
                
                
                if len(xscripts) > 1:                        
                        ujcID = row['ujc_id']
                        geneID = row['gene_id']
                        junctionStr = row['junction_string']

                        for xscript in xscripts:
                                xscriptLst.append(xscript)
                                ujcIDLst.append(ujcID)
                                geneLst.append(geneID)
                                junctionLst.append(junctionStr)
                else:
                        xscript = xscripts[0]
                        ujcID = row['ujc_id']
                        geneID = row['gene_id']
                        junctionStr = row['junction_string']

                        
                        xscriptLst.append(xscript)
                        ujcIDLst.append(ujcID)
                        geneLst.append(geneID)
                        junctionLst.append(junctionStr)


        outDf = pd.DataFrame(
                {
                        'gene_id':geneLst,
                        'transcript_id':xscriptLst,
                        'ujc_id':ujcIDLst,
                        'junction_string':junctionLst
                })
        
        return outDf

def main():
        
        """
        Run the program.

        Returns
        -------
        None.

        """
        print ("Loading...")
        omegatic = time.perf_counter()

        # if (os.path.exists(prefix + '.pickle') and os.path.getsize(prefix + '.pickle') > 0):
        #         with open(prefix + '.pickle', 'rb') as f:
        #                 exonData = pickle.load(f)
        # else:
        #         exonData = trand.io.read_exon_data_from_file(infile=args.inGTF)
        #         with open(prefix + '.pickle', 'wb') as f:
        #                 pickle.dump(exonData, f)
        
        exonData = trand.io.read_exon_data_from_file(infile=args.inGTF)
        
        toc = time.perf_counter()
        print(f"GTF Read complete! Took {toc-omegatic:0.4f} seconds. Extracting junctions...")
        tic = time.perf_counter()
        
        ujcDct = extractJunction(exonData)
                
        toc = time.perf_counter()
        print (f"Complete! Operation took {toc-tic:0.4f} seconds. Creating UJC DataFrame...")
        tic = time.perf_counter()
        
        ujcDf = createUJCDf(ujcDct=ujcDct)
        
        # if (os.path.exists(prefix + '_allUJC.pickle') and os.path.getsize(prefix + '_allUJC.pickle') > 0):
        #         with open(prefix + '_allUJC.pickle', 'rb') as f:
        #                 ujcDf = pickle.load(f)
        # else:
        #         ujcDf = createUJCDf(ujcDct=ujcDct)
        #         with open(prefix + '_allUJC.pickle', 'wb') as f:
        #                 pickle.dump(ujcDf, f)
                
        toc = time.perf_counter()
        print (f"Complete! Operation took {toc-tic:0.4f} seconds. Writing files...")
        tic = time.perf_counter()
        
        outDf = createOutput(ujcDf=ujcDf, ujcDct=ujcDct)
        
        outputPath = args.outdir + "/" + args.prefix + "_UJC_ID.csv"
       
        
        
        try:
                outDf.to_csv(outputPath, index=False)
        except OSError:
                raise OSError("Output directory must already exist.")
        
        if args.outGTF:
                
                gtfDf = createExonOutput(ujcDf=ujcDf, ujcDct=ujcDct)
                
                
                gtfOutPath = args.outdir + "/" + args.prefix + "_UJC.gtf"
                
                
                if os.path.isfile(gtfOutPath):
                        os.remove(gtfOutPath)
                
                trand.io.write_gtf(data=gtfDf, out_fhs={"gtf":gtfOutPath}, fh_name="gtf")
                

        toc = time.perf_counter()
        print(f"Complete! Operation took {toc-omegatic:0.4f} total seconds.")
        

if __name__ == '__main__':
        global args
        args = getOptions()
        
        main()
        
        

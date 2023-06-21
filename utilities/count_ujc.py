#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 13:34:55 2023

@author: k.bankole
"""

"""

Count the transcripts PER unique junction chain (UJC) of a GTF file.

Created from a previous utility in TranD named consolidation

Version 1

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
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Counts the transcripts per unique junction chain (UJC) found within a "
                                         "GTF file and outputs a summary file with this information. "
                                         "Input a GTF file (--gtf), an optional prefix for the ujc_ids (--transcript-prefix), "
                                         "an output path (--outdir) and a prefix for the output files (--prefix). "
                                         "Please note: The ujc_id is a representative transcript for all the "
                                         "transcripts with that same junction chain. The utility sorts the group of "
                                         "transcripts alphabetically and selects the first one as the representative.")
        
        ## INPUT
        parser.add_argument(
                "-g",
                "--gtf",
                dest="inGTF",
                required=True,
                help="Input a GTF file."
        )
        
        parser.add_argument(
                "-x",
                "--transcript-prefix",
                dest="trPrefix",
                required=False,
                default=None,
                help="Input a transcript prefix for the ujc_id. Useful for marking "
                        "which gtf a transcript came from. Defaults to no prefix. "
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
        
        # This takes less than 5 seconds for the small file. 
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

def createUJCDf(ujcDct, trPrefix):
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
                # it works :)
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
        
        if trPrefix:
                allUJC["ujc_id"] = trPrefix + "_" + allUJC["split"].str[0]
        else:
                allUJC["ujc_id"] = allUJC["split"].str[0]        
        allUJC.drop(columns="split")
                
        return allUJC

def createOutput(ujcDf, ignoreGene):
        """

        Parameters
        ----------
        ujcDf : DATAFRAME
                Dataframe with information on the UJCs, with their ids, transcripts, etc.
        
        ignoreGene: BOOLEAN
                User option on ignoring genes.
        
        Returns
        -------
        outDf : DATAFRAME
                A dataframe that is the list of UJCs and the number of transcripts
                within that UJC.

        """
        ujcIDLst = []
        numXscriptLst = []
        geneIDLst = []
        jStringLst = []
        
        for row in ujcDf.to_dict('records'):
                ujcID = row['ujc_id']
                xscriptStr = row['transcript_id']
                jString = row['junction_string']
                geneID = row['gene_id']
                
                numXscripts = len(xscriptStr.split('|'))
                
                ujcIDLst.append(ujcID)
                numXscriptLst.append(numXscripts)
                jStringLst.append(jString)
                geneIDLst.append(geneID)
        
        outDf = pd.DataFrame(
                {
                        'gene_id':geneIDLst,
                        'ujc_id':ujcIDLst,
                        'num_xscripts':numXscriptLst,
                        'junction_string':jStringLst
                })
        
        return outDf

if __name__ == '__main__':
        global args
        print ("Loading...")
        omegatic = time.perf_counter()
        args = getOptions()
        
        #main
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
        
        ujcDf = createUJCDf(ujcDct=ujcDct, trPrefix=args.trPrefix)
        
        # if (os.path.exists(prefix + '_allUJC.pickle') and os.path.getsize(prefix + '_allUJC.pickle') > 0):
        #         with open(prefix + '_allUJC.pickle', 'rb') as f:
        #                 ujcDf = pickle.load(f)
        # else:
        #         ujcDf = createUJCDf(ujcDct=ujcDct, trPrefix=args.trPrefix, ignoreGene=args.noGene)
        #         with open(prefix + '_allUJC.pickle', 'wb') as f:
        #                 pickle.dump(ujcDf, f)
                
        toc = time.perf_counter()
        print (f"Complete! Operation took {toc-tic:0.4f} seconds.")

        toc = time.perf_counter()
        print (f"Complete! Operation took {toc-tic:0.4f} seconds. Writing files...")
        tic = time.perf_counter()
                
        outDf = createOutput(ujcDf=ujcDf, ignoreGene=args.noGene)
        outputPath = args.outdir + "/" + args.prefix + "_UJC_count.csv"

        try:
                outDf.to_csv(outputPath, index=False)
        except OSError:
                raise OSError("Output directory must already exist.")
                
                                        
        toc = time.perf_counter()
        print(f"Complete! Operation took {toc-omegatic:0.4f} total seconds.")
        
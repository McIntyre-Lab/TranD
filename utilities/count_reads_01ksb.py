#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 13:34:55 2023

@author: k.bankole
"""

"""
(combined)
Identify the unique junction chains (UJC) of a GTF file and combine
transcripts into UJCs.

Created from a previous utility in TranD named consolidation

Version 4: Found the fastest method, separating ID and COUNT

"""

import argparse
import time
import numpy as np
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
        
        
        
        # def __str__(self):
        #         return self.seqname + ":" + str(self.lastExonEnd) + ":" + str(self.nextExonStart) + ":" + self.strand

        

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Identifies unique junction chains (UJCs) found within a "
                                         "GTF file and combines transcripts based on these UJCs. Outputs a summary file "
                                         "containing info on these combined transcripts and their \"ujc_id\". Input a GTF "
                                         "file (--infile), optional prefix for the new transcript names (--tr-prefix), "
                                         "an output path (--outdir) and a prefix for the output files (--prefix). "
                                         "Allows the option to output a new GTF file with the UJCs as the transcripts "
                                         "(--outGTF). Also allows the option ignore which gene a transcript came from"
                                         "when creating new transcript names (--ignore-gene)."
                                         "Note: \"ujc_ids\" are ranked by length. So the longest transcript under one"
                                         "UJC will prefix_gene_1, the next will be prefix_gene_2, etc. EXCEPT"
                                         "monoexons will always be number one.")
        
        ## INPUT
        parser.add_argument(
                "-i",
                "--infile",
                dest="inGTF",
                required=True,
                help="Input a GTF file."
        )
        
        parser.add_argument(
                "-t",
                "--tr-prefix",
                dest="trPrefix",
                required=False,
                default="tr",
                help="Input a prefix for the combined transcripts. Defaults to"
                "\'tr\'. (ex: tr_FBgn000000_1)"
        )
        
        parser.add_argument(
                "-ig",
                "--ignore-gene",
                dest="noGene",
                required=False,
                action="store_true",
                help="Add this argument to ignore genes when combining transcripts into UJCs."
                        "Example: \"tr_FBgn000000_1\" becomes \" tr_1.\" The gene_id in the "
                        "output GTF file will be the same as the transcript."
        )
                
        parser.add_argument(
                "-g",
                "--outGTF",
                dest="outGTF",
                action="store_false",
                help="Use this argument to remove the output of a GTF with the UJCs as transcripts. "
                        "Defaults to outputting the GTF."
        )
        
        # What does this do?
        # parser.add_argument(
        #         "-v",
        #         "--verbose",
        #         dest="verbose",
        #         action="store_true",
        #         help="dunno."
        # )
        
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
                help="Prefix for the output file(s). Example: prefix_UJC_key.csv"
        )
        
        args = parser.parse_args()
        return args

def checkStrandAndChromosome(exonData):
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
        exonDf = checkStrandAndChromosome(exonData=exonData)
        
        print ("Number of transcripts: ", end="")
        print (len(exonDf['transcript_id'].unique()))
        
        # First, instead of grouping, then sorting
        # Sort by transcript -> sort by start. the whole dataframe
        # instead of wasting time repeatedly sorting every group in the loop.
        # Which makes it now take a second. 
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
        # Takes (drum roll please...) a minute!! :)
        # Dear LORD this is fast!!!!
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
                
                # if (x == 'FBtr0077452' or x == 'FBtr0335136' or x == 'FBtr0345285'):
                #         print (info[0])
        
        return ujcDct

def createUJCDf(ujcDct, trPrefix, ignoreGene):
        monoExons = dict()
        multiExons = dict()        
        
        for xscript, info in ujcDct.items():
                junctionChain = info[7]
                
                if junctionChain:
                        newInfo = [junctionChain.chainToStr(), info[1], info[2], info[3], info[4], info[5], info[6]]
                        multiExons.update({xscript: newInfo})
                else:
                        newInfo = ['', info[1], info[2], info[3], info[4], info[5], info[6]]
                        monoExons.update({xscript: newInfo})
                        
                
        if len(monoExons) > 0:
                monoXscripts = pd.DataFrame(monoExons,
                                                index = pd.Index(["junction_string",
                                                                  "transcript_id",
                                                                  "gene_id",
                                                                  "seqname",
                                                                  "start", "end", "strand"])
                                                ).T.sort_values(by=["start", "end"])
        
                overlap = (monoXscripts['start'].shift(-1) < monoXscripts['end']).cumsum()
                
                jStringLst = []
                for row in monoXscripts.to_dict('records'):
                        jStringLst.append("monoexon_" 
                                          + str(row['start']) + "_" 
                                          + str(row['end']))
                        
                monoXscripts['junction_string'] = jStringLst
                
                monoUJC = monoXscripts.groupby(["gene_id", "junction_string"]).agg({
                        "seqname":"first",
                        "start":"min",
                        "end":"max",
                        "strand":"first",
                        "transcript_id": lambda x: '|'.join(x)}).reset_index()
                
        else:
                monoExons = None
        
        if len(multiExons) > 0:
                multiUJC = pd.DataFrame(multiExons, index=pd.Index(["junction_string", 
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
                multiExons = None
                
        if monoExons and multiExons:
                allUJC = pd.concat([monoUJC, multiUJC], ignore_index=True)
        elif monoExons:
                allUJC = monoUJC.copy()
                del(monoUJC)
        else:
                allUJC = multiUJC.copy()
                del(multiUJC)
                
        
        allUJC["ujc_length"] = allUJC["end"] - allUJC["start"]
        
        sort_order = {"gene_id": "asc", "ujc_length": "asc", "start": "asc", "transcript_id": "asc"}
        allUJC = allUJC.sort_values(by=list(sort_order.keys()), ascending=[True if val=="asc" else False for val in sort_order.values()])
        
        if not ignoreGene: 
                allUJC["transcript_rank_in_gene"] = (
                    allUJC.groupby("gene_id")["ujc_length"].rank(method="first")
                )
                
                allUJC["ujc_id"] = (
                    trPrefix
                    + "_"
                    + allUJC["gene_id"]
                    + "_"
                    + allUJC["transcript_rank_in_gene"].astype(int).map(str)
                )
        else:
                allUJC["transcript_rank_in_gene"] = (
                    allUJC["ujc_length"].rank(method="first")
                )
                
                allUJC["ujc_id"] = (
                    trPrefix
                    + allUJC['seqname']
                    + "_" 
                    + allUJC['start']
                    + "_"
                    + allUJC['end']
                    + "_"
                    + allUJC["transcript_rank_in_gene"].astype(int).map(str)
                )
                
                allUJC["gene_id"] = allUJC["ujc_id"]
                
        return allUJC

def createExonOutput(ujcDf, ujcDct):
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
        xscriptLst = []
        ujcIDLst = []
        junctionStrLst = []
        
        for row in ujcDf.to_dict('records'):
                xscripts = row['transcript_id'].split('|')
                
                
                if len(xscripts) > 1:                        
                        ujcID = row['ujc_id']
                        junctionStr = row['junction_string']
                        
                        for xscript in xscripts:
                                xscriptLst.append(xscript)
                                ujcIDLst.append(ujcID)
                                junctionStrLst.append(junctionStr)     
                        
                else:
                        xscript = xscripts[0]
                        ujcID = row['ujc_id']
                        junctionStr = row['junction_string']
                        
                        xscriptLst.append(xscript)
                        ujcIDLst.append(ujcID)
                        junctionStrLst.append(junctionStr)
                        
                
        
        outDf = pd.DataFrame(
                {
                        'transcript_id':xscriptLst,
                        'ujc_id':ujcIDLst,
                        'junctionStr':junctionStrLst
                })
        
        return outDf

if __name__ == '__main__':
        global args
        print ("Loading...")
        omegatic = time.perf_counter()
        args = getOptions()
        prefix= args.prefix
        
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
        
        ujcDf = createUJCDf(ujcDct=ujcDct, trPrefix=args.trPrefix, ignoreGene=args.noGene)
        
        # if (os.path.exists(prefix + '_allUJC.pickle') and os.path.getsize(prefix + '_allUJC.pickle') > 0):
        #         with open(prefix + '_allUJC.pickle', 'rb') as f:
        #                 ujcDf = pickle.load(f)
        # else:
        #         ujcDf = createUJCDf(ujcDct=ujcDct, trPrefix=args.trPrefix, ignoreGene=args.noGene)
        #         with open(prefix + '_allUJC.pickle', 'wb') as f:
        #                 pickle.dump(ujcDf, f)
                
        toc = time.perf_counter()
        print (f"Complete! Operation took {toc-tic:0.4f} seconds.")

        if args.outGTF:
                print ("Creating GTF file...")
                tic = time.perf_counter()
                
                gtfDf = createExonOutput(ujcDf=ujcDf, ujcDct=ujcDct)
                
                toc = time.perf_counter()
                print(f"Complete! Took {toc-tic:0.4f} seconds. Writing files...")
                
                outDf = createOutput(ujcDf=ujcDf, ujcDct=ujcDct)
                
                if not args.noGene:
                        gtfOutPath = args.outdir + "/" + args.prefix + "_UJC.gtf"
                        outputPath = args.outdir + "/" + args.prefix + "_UJC_ID.csv"
                else:
                        gtfOutPath = args.outdir + "/" + args.prefix + "_UJC_ignoregene.gtf"
                        outputPath = args.outdir + "/" + args.prefix + "_ignoregene_UJC_ID.csv"
                
                
                if os.path.isfile(gtfOutPath):
                        os.remove(gtfOutPath)
                        
                        
                trand.io.write_gtf(data=gtfDf, out_fhs={"gtf":gtfOutPath}, fh_name="gtf")
                
                outDf.to_csv(outputPath, index=False)
                
        else:
                toc = time.perf_counter()
                print(f"Complete! Took {toc-tic:0.4f} seconds. Writing files...")
                
                outDf = createOutput(ujcDf=ujcDf, ujcDct=ujcDct)
                
                if not args.noGene:
                        outputPath = args.outdir + "/" + args.prefix + "_UJC_ID.csv"
                else:
                        keyOutPath = args.outdir + "/" + args.prefix + "_ignoregene_UJC_ID.csv"
                        

                outDf.to_csv(outputPath, index=False)
                
                
        toc = time.perf_counter()
        print(f"Complete! Operation took {toc-omegatic:0.4f} total seconds.")
        
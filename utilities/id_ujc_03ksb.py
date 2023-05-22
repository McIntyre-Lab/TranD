#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 09:31:01 2023

@author: k.bankole
"""

"""

Identify the unique junction chains (UJC) of a GTF file and combine
transcripts into UJCs.

Created from a previous utility in TranD named consolidation

Version 3: Playing with generators!

"""

import argparse
import time
import pandas as pd
import trand.io
import os
from dataclasses import dataclass
import pickle
import copy
from itertools import islice
import numpy as np

# TODO list:
        # verbose?
        # option to not ouptut the GTF
        
@dataclass()
class JUNCTION:
        seqname: str
        lastExonEnd: int
        nextExonStart: int
        strand: str
        
        def __str__(self):
                return self.seqname + ":" + str(self.lastExonEnd) + ":" + str(self.nextExonStart) + ":" + self.strand

@dataclass()
class J_CHAIN:
        
        junctionLst: list
        
        def chainToStr(self):
                chainStr = []
                
                for junction in self.junctionLst:
                        chainStr.append(str(junction))
                
                return '|'.join(chainStr)
        
        def __hash__(self):
                return hash(self.chainToStr())
        
        def __eq__(self, other):
                return self.junctionLst == other.junctionLst
        
        
def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Identifies unique junction chains (UJCs) found within a "
                                         "GTF file and combines transcripts based on these UJCs. Outputs a summary key file "
                                         "Containing info on these combined transcripts. Input a GTF "
                                         "file (--infile), optional prefix for the new transcript names (--tr-prefix), "
                                         "an output path (--outdir) and a prefix for the output files (--prefix). "
                                         "Allows the option to output a new GTF file with the UJCs as the transcripts "
                                         "(--outGTF). Also allows the option ignore which gene a transcript came from"
                                         "when creating new transcript names (--ignore-gene)")
        
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

# def generateObject(df):
        
#         for row in df.values.tolist()[:-1]:
#                 junction = JUNCTION(seqname=row[0], lastExonEnd=str(row[2]), 
#                                     nextExonStart=row[1], strand=row[3])
                
#                 yield junction
        
        
        

# this seems to be the fastest way. list > .values > tuples > todict
def extractJunction(exonData):
        exonData = checkStrandAndChromosome(exonData=exonData)
        xscriptGrps = exonData.groupby("transcript_id")
        
        ujcDct = {}
        
        for xscriptID, group in xscriptGrps:
                
                junctionLst = []
                if len(group) > 1:
                        sortedGrp = group.sort_values(by="start").reset_index(drop=True)
                        
                        sortedGrp['start'] = sortedGrp['start'].shift(-1).fillna(0).astype(int).astype(str)
                        
                        columns = sortedGrp['seqname'], sortedGrp['end'], sortedGrp['start'], sortedGrp['strand']
                        
                        junctionLst = [JUNCTION(seqname=seqname, lastExonEnd=end, nextExonStart=start, strand=strand) 
                                        for seqname, end, start, strand in islice(zip(*columns), len(sortedGrp) - 1)]
                        
                        start = group['start'].min()
                        end = group['end'].max()
                        
                        junctionChain = J_CHAIN(junctionLst)
                else:
                        junctionChain = None
                        start = group['start'].iloc[0]
                        end = group['end'].iloc[0]
                
                                        
                seqname = group['seqname'].iloc[0]
                strand = group['strand'].iloc[0]
                geneID = group['gene_id'].iloc[0]
                
                ujcDct[xscriptID] = [junctionChain,
                                     xscriptID,
                                     geneID,
                                     seqname,
                                     start,
                                     end,
                                     strand]
                
        return ujcDct
                
def createUJCDf(ujcDct, trPrefix, ignoreGene):
        print (len(ujcDct))
        monoExons = {t:ujcDct[t] for t in ujcDct if not ujcDct[t][0]}
        if len (monoExons) > 0:
                monoXscripts = pd.DataFrame(monoExons,
                                                index = pd.Index(["junction_string",
                                                                 "transcript_id",
                                                                 "gene_id",
                                                                 "seqname",
                                                                 "start", "end", "strand"])
                                                ).T.sort_values(by=["start", "end"])
                
                overlap = (monoXscripts['start'].shift(-1) < monoXscripts['end']).cumsum()
                monoXscripts['junction_string'] = overlap + 1
                
                monoUJC = monoXscripts.groupby(["gene_id", "junction_string"]).agg({
                        "seqname":"first",
                        "start":"min",
                        "end":"max",
                        "strand":"first",
                        "transcript_id": lambda x: '|'.join(x)}).reset_index()
                
                del (monoXscripts)
        else:
                monoExons = None
        
        
        multiExons = copy.deepcopy(ujcDct)
        multiExons = {t:multiExons[t] for t in multiExons if multiExons[t][0]}

        for value in multiExons.values():
                value[0] = value[0].chainToStr()
        

        # del (ujcDct?)
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
        
        for row in ujcDf.to_dict('records'):
                # only takes the first transcript id in the list if there are multiple
                # this should be fine? if the transcript_id is listed under the same junction chain
                # they have the same junctions -> same exons?
                                
                xscript = row['transcript_id'].split('|')[0]
                seqname = row['seqname']
                start = row['start']
                end = row['end']
                strand = row['strand']
                ujcID = row['ujc_id']
                geneID = row['gene_id']
                

                # if it is not monoexon (has a junctionLst)
                if (ujcDct.get(xscript)[0]):
                        junctions = ujcDct.get(xscript)[0].junctionLst
                                                
                        # donors = end of last exon 
                        # accceptor = start of next exon
                        startLst.append(start)
                        endLst.append(junctions[0].lastExonEnd)
                        seqnameLst.append(seqname)
                        ujcIDLst.append(ujcID)
                        geneIDLst.append(geneID)
                        strandLst.append(strand)
                        
                        #all the ends are shifted down 1
                        
                        for idx in range(0,len(junctions)-1):
                                endLst.append(junctions[idx+1].lastExonEnd)
                                startLst.append(junctions[idx].nextExonStart)
                                seqnameLst.append(seqname)
                                ujcIDLst.append(ujcID)
                                geneIDLst.append(geneID)
                                strandLst.append(strand)
                        
                        startLst.append(junctions[-1].nextExonStart)
                        endLst.append(end)
                        seqnameLst.append(seqname)
                        ujcIDLst.append(ujcID)
                        geneIDLst.append(geneID)
                        strandLst.append(strand)
                
                else:
                        seqnameLst.append(seqname)
                        startLst.append(start)
                        endLst.append(end)
                        strandLst.append(strand)
                        ujcIDLst.append(ujcID)
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

def splitColBySep(df, colName=None, sep=None, sortLst=None):
        # Split a column of a DF by some charaacter like '|' and keep all other values the same
        
        if colName == None:
                colName = 'transcript_id'
        
        if sep == None:
                sep = '|'
        
        splitLst = df[colName].str.split(sep).apply(pd.Series, 1).stack()
        splitLst.index = splitLst.index.droplevel(-1)
        
        tempDf = df.copy()
        
        del (tempDf[colName])
        splitDf = tempDf.join(splitLst.rename(colName))
        
        if sortLst != None:
                splitDf = splitDf.sort_values(by=sortLst)
        
        del (tempDf, splitLst)
        return splitDf


def createKeyFile(ujcDf):
        keyDf = splitColBySep(df=ujcDf[['gene_id', 'transcript_id', 'ujc_id']],colName='transcript_id', sep = '|')
        
        keyDf['num_xscripts_in_ujc'] = keyDf.groupby(
                'ujc_id')['transcript_id'].transform('nunique')
        
        keyDf['num_xscript_in_gene'] = keyDf.groupby(
                'gene_id')['transcript_id'].transform('nunique')
        
        keyDf['num_ujc_in_gene'] = keyDf.groupby(
                'gene_id')['ujc_id'].transform('nunique')
        
        keyDf['flag_gene_consolidated'] = np.where(
                keyDf['num_xscript_in_gene'] > keyDf["num_ujc_in_gene"],
                1,
                0)
        
        
        return keyDf
        

# move everything here eventually
def main():
        return 'KSB'

if __name__ == '__main__':
        global args
        
        print ("Loading...")
        omegatic = time.perf_counter()
        args = getOptions()
        
        # main
        if (os.path.exists('exonData.pickle') and os.path.getsize('exonData.pickle') > 0):
                with open('exonData.pickle', 'rb') as f:
                        exonData = pickle.load(f)
        else:
                exonData = trand.io.read_exon_data_from_file(infile=args.inGTF)
                with open('exonData.pickle', 'wb') as f:
                        pickle.dump(exonData, f)        
        
        
        toc = time.perf_counter()
        print(f"GTF Read complete! Took {toc-omegatic:0.4f} seconds. Extracting junctions...")
        
        tic = time.perf_counter()
                        
        # if (os.path.exists('ujcDct.pickle') and os.path.getsize('ujcDct.pickle') > 0):
        #         with open('ujcDct.pickle', 'rb') as f:
        #                 masterUJCDct = pickle.load(f)
        # else:
        masterUJCDct = extractJunctiona(exonData=exonData)
        #         with open('ujcDct.pickle', 'wb') as f:
        #                 pickle.dump(masterUJCDct, f)  

        toc = time.perf_counter()
        print(f"Junction extraction complete! Took {toc-tic:0.4f} seconds. Creating UJC Df...")
        
        tic = time.perf_counter()

        ujcDf = createUJCDf(ujcDct=masterUJCDct, trPrefix=args.trPrefix, ignoreGene=args.noGene)
        

        
        if args.outGTF:
                toc = time.perf_counter()
                print(f"Complete! Took {toc-tic:0.4f} seconds. Creating exon output...")
                
                tic = time.perf_counter()
                
                exonDf = createExonOutput(ujcDf=ujcDf, ujcDct=masterUJCDct)
                
                toc = time.perf_counter()
                print(f"Complete! Took {toc-tic:0.4f} seconds. Writing files...")
                
                if not args.noGene:
                        gtfOutPath = args.outdir + "/" + args.prefix + "_UJC.gtf"
                        keyOutPath = args.outdir + "/" + args.prefix + "_UJC_key.csv"
                else:
                        gtfOutPath = args.outdir + "/" + args.prefix + "_UJC_ignoregene.gtf"
                        keyOutPath = args.outdir + "/" + args.prefix + "_UJC_ignoregene_key.csv"
                
                
                if os.path.isfile(gtfOutPath):
                        os.remove(gtfOutPath)
                        
                        
                trand.io.write_gtf(data=exonDf, out_fhs={"gtf":gtfOutPath}, fh_name="gtf")
                
                
                keyDf = createKeyFile(ujcDf=ujcDf)
                keyDf.to_csv(keyOutPath, index=False)
                
        else:
                toc = time.perf_counter()
                print(f"Complete! Took {toc-tic:0.4f} seconds. Writing files...")
                
                if not args.noGene:
                        keyOutPath = args.outdir + "/" + args.prefix + "_UJC_key.csv"
                else:
                        keyOutPath = args.outdir + "/" + args.prefix + "_UJC_ignoregene_key.csv"
                
                keyDf = createKeyFile(ujcDf=ujcDf)
                keyDf.to_csv(keyOutPath, index=False)
                
                
        toc = time.perf_counter()
        print(f"Complete! Operation took {toc-omegatic:0.4f} total seconds.")
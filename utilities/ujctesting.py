#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 12:07:12 2023

@author: k.bankole
"""

import argparse
import time
import trand.io
import pandas as pd
from operator import itemgetter
import numpy as np
from dataclasses import dataclass
import copy

import pickle
import os

# MISSING STUFF: verbose, actual argparse arguments
# Potential Ideas that I need a whole day for:
        # reading the GTF directly instead of usign read_exon (extract the info
        # I need into a DF instead of just using what TranD uses)
        
        

@dataclass()
class JUNCTION:
        seqname: str
        lastExonEnd: int
        nextExonStart: int
        strand: str
        
        def __str__(self):
                return (self.seqname + ":" + str(self.lastExonEnd) + ":" + str(self.nextExonStart) + ":" + str(self.strand))
        
@dataclass()
class J_CHAIN:
        # Enter a sorted junction list, the following are a bunch of transformation,
        # comparison, and string representation functions
        
        junctionLst: list
        
        def chainToStr(self):
                
                chainStr = [];
                
                for junction in self.junctionLst:
                        chainStr.append(str(junction))
                
                return '|'.join(chainStr)
        
        
        

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Consolidate GTF by identifying unique junction chains.")
        
        # Input data
        # parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="Input GTF file.")
        # parser.add_argument("-c", "--consol-prefix", dest="consolPrefix", required=False, default="tr", help="Prefix to use for consolidation transcript id values.")
        # parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_true", help="Verbose mode.")
        
        # # Output data
        # parser.add_argument("-o", "--output-prefix", dest="outPrefix", required=True, help="Full path and prefix for all output files.")
        
        parser.add_argument(
                "-c",
                "--choose",
                dest="yeah",
                required=True,
                help="No."
        )
        
        args = parser.parse_args()
        return args

# def main():
#         return None

# all good, takes <10s on mel2mel (3million exons)
def check_strand_and_chromosome(exon_data):
        gene_groups = exon_data.groupby("gene_id")
        strand_check = gene_groups["strand"].nunique()
        
        if (strand_check > 1).any():
            bad_strand = list(strand_check[strand_check>1].index)
            for gene in bad_strand:
                print("!!! WARNING: gene {} contains transcripts/exons on both strands - removing from consolidation".format(gene))
            exon_data = exon_data[~exon_data["gene_id"].isin(bad_strand)]
        chr_check = gene_groups["seqname"].nunique()
        if (chr_check > 1).any():
            bad_chr = list(chr_check[chr_check>1].index)
            for gene in bad_chr:
                print("!!! WARNING: gene {} contains transcripts/exons on difference chromosomes - removing from consolidation".format(gene))
            exon_data = exon_data[~exon_data["gene_id"].isin(bad_chr)]
        return exon_data

# for time comp
def oldExtractJunction(exon_data):
        """Create a junction strings for each transcript"""
        
        # Check for genes with inconsistent chromsome or strand
        exon_data = check_strand_and_chromosome(exon_data)
        
        transcript_groups = exon_data.groupby("transcript_id")
        
        # Iterate over transcript groups
        junction_chains = {}
        # num = 0
        for transcript_id, group in transcript_groups:
                # num = num + 1
                # if num > 100:
                #     break
                if len(group) > 1:
                        # Sort exons by start position
                        sorted_exons = group.sort_values(by="start").reset_index(drop=True)
                        
                        # Create junction strings using vectorized operations
                        junctions = sorted_exons["seqname"].shift(-1).str.cat(
                                sorted_exons["end"].astype(int).astype(str), sep=":", na_rep=""
                            ).str.cat(
                                sorted_exons["start"].shift(-1).fillna(0).astype(int).astype(str), sep=":", na_rep=""
                            ).str.cat(
                                sorted_exons["strand"], sep=":", na_rep=""
                            ).tolist()[:-1]
                        
                        # Concatenate junctions into a single string
                        junction_chain = '|'.join(junctions)
                        start = sorted_exons["start"].min()
                        end = sorted_exons["end"].max()
                else:
                        junction_chain = ""
                        start = group["start"].iat[0]
                        end = group["end"].iat[0]
            
                strand = group["strand"].iat[0]
                seqname = group["seqname"].iat[0]
                gene_id = group["gene_id"].iat[0]
                junction_chains[transcript_id] = [junction_chain, 
                                                  transcript_id, 
                                                  gene_id, 
                                                  seqname, 
                                                  start, 
                                                  end, 
                                                  strand]
       
        return junction_chains

def extractJunction(exon_data): 
        exon_data = check_strand_and_chromosome(exon_data=exon_data)
        
        transcript_groups = exon_data.groupby("transcript_id")
        
        # note to self: a UJC is just a transcript, a collection of junctions
        
        ujcDct = {}
        
        for transcript_id, group in transcript_groups:
                
                junctionLst = []
                if len(group) > 1:
                        sortedgroup = group.sort_values(by='start').reset_index(drop=True)
                        
                        # Shifts the entire column up one -> next exon start
                        sortedgroup['start'] = sortedgroup['start'].shift(-1).fillna(0).astype(int).astype(str)
                        
                        #print (sortedgroup)                        
                        for row in sortedgroup.to_dict('records')[:-1]:
                                seqname = row['seqname']
                                strand = row['strand']
                                nextExonStart = row['start']
                                lastExonEnd = str(row['end'])
                        
                                junctionLst.append(JUNCTION(seqname=seqname, lastExonEnd=lastExonEnd, 
                                                            nextExonStart=nextExonStart, strand=strand))
                                
                        start = group['start'].min()
                        end = group['end'].max()
                                                                
                        junctionChain = J_CHAIN(junctionLst)
                        
                else:
                        junctionChain = None
                        start = group["start"].iat[0]
                        end = group["end"].iat[0]
                
                seqname = group["seqname"].iat[0]
                strand = group["strand"].iat[0]
                gene_id = group["gene_id"].iat[0]
                
                ujcDct[transcript_id] = [junctionChain,
                                         transcript_id, 
                                         gene_id, seqname, start, end, strand]
                
                
        return ujcDct

# change name
def createAllUJC(ujcDct):
        
        monoExons = {t:ujcDct[t] for t in ujcDct if not ujcDct[t][0]}
        # monoExons = {t:oldUjcDct[t] for t in oldUjcDct if oldUjcDct[t][0] == ""}
        if len(monoExons) > 0:
                monoexon_transcripts = pd.DataFrame(monoExons, 
                                                    index=pd.Index(["junction_string", 
                                                                    "transcript_id", 
                                                                    "gene_id", 
                                                                    "seqname", 
                                                                    "start", "end", "strand"])
                                                    ).T.sort_values(by=["start", "end"])
                
                # compares the start of one to the end of the other for overlap -> counts
                # each junction chain is just a number
                overlap = (monoexon_transcripts["start"].shift(-1) < monoexon_transcripts["end"]).cumsum()
                monoexon_transcripts["junction_string"] = overlap + 1
                
                
                monoUJC = monoexon_transcripts.groupby(["gene_id","junction_string"]).agg({
                    "seqname": "first",
                    "start": "min",
                    "end": "max",
                    "strand": "first",
                    "transcript_id": lambda x: "|".join(x)}).reset_index()
                
                del (monoexon_transcripts)
        else:
                monoExons = None
        
        
        multiExons = copy.deepcopy(ujcDct)
        multiExons = {t:multiExons[t] for t in multiExons if multiExons[t][0]}
        # oldmultiExons = {t:oldUjcDct[t] for t in oldUjcDct if oldUjcDct[t][0] != ""}

        # del (ujcDct)
        
        for value in multiExons.values():
                value[0] = value[0].chainToStr()
        
        if len(multiExons) > 0:
                multiUJC = pd.DataFrame(multiExons, index=pd.Index(["junction_string", 
                                                                    "transcript_id", 
                                                                    "gene_id", 
                                                                    "seqname", "start", "end", "strand"])
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
                del(monoUJC)
        
        allUJC["ujc_length"] = allUJC["end"] - allUJC["start"]
        
        
        sort_order = {"gene_id": "asc", "ujc_length": "asc", "start": "asc", "transcript_id": "asc"}
        allUJC = allUJC.sort_values(by=list(sort_order.keys()), ascending=[True if val=="asc" else False for val in sort_order.values()])
        
        allUJC["transcript_rank_in_gene"] = (
            allUJC.groupby("gene_id")["ujc_length"].rank(method="first")
        )
        
        allUJC["ujc_id"] = (
            "tr"
            + "_"
            + allUJC["gene_id"]
            + "_"
            + allUJC["transcript_rank_in_gene"].astype(int).map(str)
        )
        
        return allUJC
        

# WORKS!!!!!!
def createExonOutput(df, ujcDct):
        
        seqnameLst = []
        startLst = []
        endLst = []
        strandLst = []
        ujcIDLst = []
        geneIDLst = []
        
        for row in df.to_dict('records'):
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

def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
        # Split variable by some character like '|' or ',' and keep all other values the same
        if col_name == None:
                col_name = 'transcript_id'
        if sep == None:
                sep = "|"
                
        
        splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
        splitList.index = splitList.index.droplevel(-1)
        
        tempDF = df.copy()
        del(tempDF[col_name])
        
        splitDF = tempDF.join(splitList.rename(col_name))
        if sort_list != None:
                splitDF = splitDF.sort_values(by=sort_list)
        del(tempDF, splitList)
        
        return splitDF


if __name__ == '__main__':
        # Parse command line arguments
        global args
        
        tic = time.perf_counter()
        args = getOptions()
        
        path = ("/ufgi-vulcana/data/k.bankole/github/TranD/testdata/dmelexondata.csv" if args.yeah.upper() == "Y"
                else "/ufgi-vulcana/data/k.bankole/github/TranD/testdata/mel2melexondata.csv")
                
        
        exon_data = pd.read_csv(path, low_memory=False)
                
        if (os.path.exists('ujcDct.pickle') and os.path.getsize('ujcDct.pickle') > 0):
                with open('ujcDct.pickle', 'rb') as f:
                        ujcDct = pickle.load(f)
        else:
                ujcDct = extractJunction(exon_data=exon_data)
                with open('ujcDct.pickle', 'wb') as f:
                        pickle.dump(ujcDct, f)
        
        if (os.path.exists('oldUjcDct.pickle') and os.path.getsize('oldUjcDct.pickle') > 0):
                with open('oldUjcDct.pickle', 'rb') as f:
                        oldUjcDct = pickle.load(f)
                        
        else:
                oldUjcDct = oldExtractJunction(exon_data=exon_data)
                with open('oldUjcDct.pickle', 'wb') as f:
                        pickle.dump(oldUjcDct, f)
                        

        
        if (os.path.exists('allUJC.pickle') and os.path.getsize('allUJC.pickle') > 0):
                with open('allUJC.pickle', 'rb') as f:
                        allUJC = pickle.load(f)
                        
        else:
                allUJC = createAllUJC(ujcDct)
                with open('allUJC.pickle', 'wb') as f:
                        pickle.dump(allUJC, f)
        
        exonDf = createExonOutput(df=allUJC, ujcDct=ujcDct)        
        
        keys = split_column_by_sep(allUJC[["gene_id", "transcript_id", "ujc_id"]], col_name="transcript_id", sep="|")
        
        
        keys["num_transcript_in_consol_transcript"] = keys.groupby(
                "ujc_id")["transcript_id"].transform('nunique')
        
        keys["num_transcript_id_in_gene"] = keys.groupby(
                "gene_id")["transcript_id"].transform('nunique')
        
        keys["num_consol_transcript_id_in_gene"] = keys.groupby(
                "gene_id")["ujc_id"].transform('nunique')
        
        keys["flag_gene_consolidated"] = np.where(
            keys["num_transcript_id_in_gene"] > keys["num_consol_transcript_id_in_gene"],
            1,
            0
        )
        keys.to_csv("testkey.csv")
        
        trand.io.write_gtf(exonDf, {"gtf":"outputgtf.gtf"}, "gtf")
        
        toc = time.perf_counter()
        print(f"Complete! Operation took {toc-tic:0.4f} seconds.")
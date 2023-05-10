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

import pickle
import os

# MISSING STUFF: verbose, actual argparse arguments
# Potential Ideas that I need a whole day for:
        # reading the GTF directly instead of usign read_exon (extract the info
        # I need into a DF instead of just using what TranD uses)
        
        

@dataclass()
class JUNCTION:
        seqname: str
        start: int
        end: int
        strand: str
        
        def __str__(self):
                return (self.seqname + ":" + str(self.start) + ":" + str(self.end) + ":" + str(self.strand))
        
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
                        
                        for row in group.to_dict('records'):
                                seqname = row['seqname']
                                strand = row['strand']
                                jcStart = row['start']
                                jcEnd = row['end']
                                
                                junctionLst.append(JUNCTION(seqname, jcStart, jcEnd, strand))
                                
                        start = group['start'].min()
                        end = group['end'].max()
                
                        junctionLst.sort(key=lambda x: x.start)
                        
                        #junctionChain = J_CHAIN(junctionLst)
                        #break;
                else:
                        junctionLst = None
                        start = group["start"].iat[0]
                        end = group["end"].iat[0]
                
                seqname = group["seqname"].iat[0]
                strand = group["strand"].iat[0]
                gene_id = group["gene_id"].iat[0]
                
                ujcDct[transcript_id] = [junctionLst,
                                         transcript_id, 
                                         gene_id, seqname, start, end, strand]
        return ujcDct

        

if __name__ == '__main__':
        # Parse command line arguments
        global args
        global ujcDct
        
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
                        

        monoExons = {t:ujcDct[t] for t in ujcDct if not ujcDct[t][0]}
        monoExons = {t:oldUjcDct[t] for t in oldUjcDct if oldUjcDct[t][0] == ""}
        if len(monoExons) > 0:
                oldmonoexon_transcripts = pd.DataFrame(monoExons, 
                                                    index=pd.Index(["junction_string", 
                                                                    "transcript_id", 
                                                                    "gene_id", 
                                                                    "seqname", 
                                                                    "start", "end", "strand"])
                                                    ).T.sort_values(by=["start", "end"])
                
                # compares the start of one to the end of the other for overlap -> counts
                # each junction chain is just a number
                overlap = (oldmonoexon_transcripts["start"].shift(-1) < oldmonoexon_transcripts["end"]).cumsum()
                oldmonoexon_transcripts["junction_string"] = overlap + 1
                
        #         monoUJC = monoexon_transcripts.groupby(["gene_id","junction_string"]).agg({
        #             "seqname": "first",
        #             "start": "min",
        #             "end": "max",
        #             "strand": "first",
        #             "transcript_id": lambda x: "|".join(x)}).reset_index()

        #         del(monoexon_transcripts)
                
        # multiExons = {t:ujcDct[t] for t in ujcDct if ujcDct[t][0]}#!= ""}

                
        
        # if len(multiExons) > 0:
        #         multiUJC = pd.DataFrame(multiExons, index=pd.Index(["junction_string", 
        #                                                             "transcript_id", 
        #                                                             "gene_id", 
        #                                                             "seqname", "start", "end", "strand"])
        #             ).T.groupby(["gene_id","junction_string"]).agg({
        #         "seqname": "first",
        #         "start": "min",
        #         "end": "max",
        #         "strand": "first",
        #         "transcript_id": lambda x: "|".join(x)}).reset_index()

        
                
        toc = time.perf_counter()
        print(f"Complete! Operation took {toc-tic:0.4f} seconds.")
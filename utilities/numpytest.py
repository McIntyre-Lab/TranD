#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 13:34:55 2023

@author: k.bankole
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
        
        # def __str__(self):
        #         return self.seqname + ":" + str(self.lastExonEnd) + ":" + str(self.nextExonStart) + ":" + self.strand

        

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
        
        df = checkStrandAndChromosome(exonData=exonData)
        
        print ("Number of transripts: ", end="")
        print (len(df['transcript_id'].unique()))
        
        # First, insted of grouping, then sorting
        # Sort by transcript -> sort by start. the whole dataframe
        # instead of wasting time repeatedly sorting every group in the loop.
        # Which makes it now take a second. lol.
        sortedDf = exonData.sort_values(by=['transcript_id', 'start']).reset_index(drop=True)
        
        xscriptGrps = exonData.groupby("transcript_id")
        
        ujcDct = {}
        
        addedXscripts = set()
        
        for row in sortedDf.to_dict('records'):                
                xscript = row['transcript_id']
                
                seqname = row['seqname']
                strand = row['strand']
                geneID = row['gene_id']
                
                start = row['start']
                end = row['end']
                
                if xscript in ujcDct:
                        exonLst = ujcDct[xscript][0]
                        
                        junctionLst.append((end,start))
                else:
                        junctionLst = [(end,start)]
                        ujcDct[xscript] = [junctionLst,
                                           xscript,
                                           geneID,
                                           seqname,
                                           start,
                                           end,
                                           strand]
                
                
        # for xscriptID, group in xscriptGrps:
        #         group = group
        
        # for xscriptID, group in xscriptGrps:
                
        #         seqname = group['seqname'].iloc[0]
        #         strand = group['strand'].iloc[0]
        #         geneID = group['gene_id'].iloc[0]
                

                
        #         if len(group) > 1:
                        
        #                 startValues = group['start'].to_numpy()
        #                 endValues = group['end'].to_numpy()
                        
        #                 # dear god this is fast. (24 seconds for the small file.)
        #                 # Drum Roll Please....... 400s for the 3 milly file
                        
        #                 junctions = list (zip(endValues[:-1], startValues[1:]))
                        
        #                 del (startValues)
        #                 del (endValues)
                        
        #                 junctionChain = J_CHAIN(seqname=seqname, junctionLst=junctions, strand=strand)
                        
        #                 start = group['start'].min()
        #                 end = group['end'].max()
        #         else:
        #                 junctionChain = None
        #                 start = group['start'].iloc[0]
        #                 end = group['end'].iloc[0] 
                
                
                
        #         ujcDct[xscriptID] = [junctionChain,
        #                               xscriptID,
        #                               geneID,
        #                               seqname,
        #                               start,
        #                               end,
        #                               strand]
                

                
        
                
                
                
                
        
        
        toc = time.perf_counter()
        print (f"Complete! Operation took {toc-tic:0.4f} seconds.")
        
        
        
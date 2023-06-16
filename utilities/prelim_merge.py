#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:04:12 2023

@author: k.bankole
"""

import argparse
import pandas as pd

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(
            description=(
                "Create a Transcript Map File from TranD output files."
            )
        )
        
        parser.add_argument(
            "-e",
            "--ERG-file",
            dest="ERG",
            required=True,
            help=(
                "Input file location for xscript_output ERG file created from TranD output data."
            )
        )
        
        parser.add_argument(
            "-p",
            "--pairwise",
            dest="PD",
            required=True,
            help=(
                "Input file location of pairwise distance file."
            )
        )
        
        parser.add_argument(
            "-1",
            "--name1",
            dest="GTF1",
            required=False,
            default="gtf1",
            help=(
                "Optional name for the first GTF when processing 2 GTF output."
                "Added to the header of the column when merged."
                "Default: gtf1"
            )
        )
        
        parser.add_argument(
            "-2",
            "--name2",
            dest="GTF2",
            required=False,
            default="gtf2",
            help=(
                "Optional name for the first GTF when processing 2 GTF output."
                "Added to the header of the column when merged."
                "Default: gtf2"
            )
        )
                
        parser.add_argument(
            "-o",
            "--outdir",
            dest="outDir",
            required=False,
            help=(
                    "Output directory, must already exist."
            )
        )
        
        args = parser.parse_args()
        return args


#def main():
        

if __name__ == '__main__':
        # Parse command line arguments
        global args
        args = getOptions()
        
        ergDf = pd.read_csv(args.ERG, low_memory=False)
        
        print ("Num ERG Xscripts: " + str(len(ergDf)))
        
        ergDf1 = ergDf[ergDf['which_gtf'] == 1]
        print ("Num GTF1 XSCRIPTS: " + str(len(ergDf1)))
        
        ergDf2 = ergDf[ergDf['which_gtf'] == 2]
        print ("Num GTF1 XSCRIPTS: " + str(len(ergDf2)))    

        pdDf = pd.read_csv(args.PD, low_memory=False)
        
        unqXscriptSet = set(pd.concat([pdDf['transcript_1'], pdDf['transcript_2']]))
        
        print ("Num PD Xscripts: " + str(len(unqXscriptSet)))        
        
        for col in pdDf.columns:
                
                if "flag_min_match_" in col:
                        gtf1 = col[len(("flag_min_match_")):]
                        break
        
        for col in pdDf.columns:
                if "flag_min_match_" in col:
                        gtf2 = col[len(("flag_min_match_")):]
                        
                        if gtf2 == gtf1:
                                continue
                        else:
                                break
                
        subsetDf = pdDf[(pdDf['flag_min_match_' + gtf1] == 1) | (pdDf['flag_min_match_' + gtf2] == 1)]
        
        subsetDf = subsetDf[subsetDf['flag_RMP'] == 1]
                
        print ("Num true minimums: " + str(len(subsetDf)))
                

    
#   main()
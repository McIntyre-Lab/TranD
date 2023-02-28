#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:43:44 2023

@author: k.bankole
"""

"""
Identify possible ER sets using TRAND ouptput of a 1 or 2 GTF pairwise file 

"""

import argparse
import pandas as pd

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        parser = argparse.ArgumentParser(description="Output a csv containing information on "
                                         "the ER sets for a set of transcripts using TRAND output "
                                         "data. Contains the option to include or exclude "
                                         "transcripts with intron retention events (--includeIR). "
                                         "Input a TRAND output file (csv) (--indir), "
                                         "IR inclusion (--includeIR), and an output path (--outdir)."
                                         )
        
        # INPUT
        parser.add_argument(
                "-i",
                "--indir",
                dest="indir",
                required=True,
                help="Location of TRAND output file")
        
        parser.add_argument(
                "-ir",
                "--includeIR",
                dest="includeIR",
                required=True,
                default='Y',
                const='Y',
                nargs='?',
                choices=['Y','N'],
                help="Choose N to exclude transcript models with intron retention events from ER sets.")
        
        # OUTPUT
        parser.add_argument(
                "-o",
                "--outdir",
                dest="outdir",
                required=True,
                help="Location of output directory, created if missing")
        
        parser.add_argument(
                "-f",
                "--force",
                dest="force",
                action="store_true",
                help="Force overwrite existing output directory and files within")
        
        args = parser.parse_args()
        return args

# input dataframe and chop. (input is a full csv file)
# necessary columns: t1, t2, flag_IR, prop_ER
# -> new dataframe: t1, t2, flag_ER (with all IR removed if excludeIR = true)




def identifyERSet(inDf, intronRetention):
        
        pairERInfo = inDf[
                [
                        "gene_id",
                        "transcript_1",
                        "transcript_2",
                        "prop_ER_similar",
                        "flag_IR"
                ]
        ].copy()

        if (not intronRetention):
                invlvdInIR = (pd.concat(pairERInfo[pairERInfo["flag_IR"==1]["transcript_1"]], 
                                                   pairERInfo[pairERInfo["flag_IR"==1]["transcript_2"]]).unique())
                exRegInfo = pairERInfo[(pairERInfo["flag_IR"] == 0)]
        # delete this post-testing
        else:
                invlvdInIR = pd.DataFrame({"boing": [1,2]})
                
        outDf = invlvdInIR
        return outDf, pairERInfo

def main():
        
        # input csv to dataframe
        inputDf = pd.read_csv(args.indir)
        
        return 'KSB'

if __name__ == '__main__':
        global args
        args = getOptions()
        
        inputDf = pd.read_csv(args.indir)

        if (args.includeIR == 'Y'):
                outputDf,pairER = identifyERSet(inputDf, True)
        elif (args.includeIR == 'N'):
                outputDf,pairER = identifyERSet(inputDf, False)
                
        
        output = pd.DataFrame.to_csv(outputDf)
        main()
        
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
        
        # unfiltered input dataframe, chopped to neccessary info
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
                # working, stores all IR flaggged xscripts in a bin ->
                invlvdInIR = pd.concat([pairERInfo[pairERInfo['flag_IR']==1]['transcript_1'],
                                       pairERInfo[pairERInfo['flag_IR']==1]['transcript_2']]).unique()
                                
                
                # working, -> remove all IR xscripts from working dataframe
                exRegInfo = pairERInfo[pairERInfo["flag_IR"] == 0]
        else: # just here for testing purposes, will remove
                print("IR Retention included")



        exRegInfo = pairERInfo[
                [
                        'transcript_1',
                        'transcript_2',
                        'prop_ER_similar'
                ]
        ]
        
        
        # ER similar proportion to boolean
        pd.options.mode.chained_assignment = None
        exRegInfo['prop_ER_similar'] = exRegInfo['prop_ER_similar'] == 1
        pd.options.mode.chained_assignment = "warn"
        
        exRegInfo = exRegInfo.rename(columns={"prop_ER_similar":"flag_ER_full_overlap"})
        
        # this could end up being useful??
        xscriptMstrLst = pd.concat([exRegInfo['transcript_1'], exRegInfo['transcript_2']]).unique()
        
        lftvrXscriptLst = xscriptMstrLst.copy().tolist()
        removedXscriptLst = []        
        # ok so iterating through pandas = terrible
        
        # iterate through a dictionary version (should be faster) 
        # structure: each row is a dictionary
        # key = column header (same for every row)
        # value = specific row value
        
        # access flag_ER: row['flag_ER_full_overlap']
        # access xscript: row['transcript_1'] or row['transcript_2']
        iterDict = exRegInfo.to_dict('records')

        # contains all exon region sets        
        exRegSetLst = []      

        # counting the number of exon region sets
        for row in iterDict:
                fullOvlpFlag = row['flag_ER_full_overlap']
                xscript1 = row['transcript_1']
                xscript2 = row['transcript_2']
                
                # should work... scripts w no overlap placed in "leftovers"
                if (fullOvlpFlag):
                        
                        #if its the first set (master list is empty)
                        if not exRegSetLst:
                                exRegSetLst.append([xscript1, xscript2])
                        else:
                                newLst = iterateAllERSet(xscript1, xscript2, exRegSetLst)
                                
                                if newLst is not None:
                                        exRegSetLst.append(newLst)
                                        
                        # removes transcript pair from leftover list (new variable name!!)
                        # removes if not already removed
                        for lst in exRegSetLst:
                                if (xscript1 in lst) and (xscript1 not in removedXscriptLst):
                                        lftvrXscriptLst.remove(xscript1)
                                        removedXscriptLst.append(xscript1)
                                if (xscript2 in lst) and (xscript2 not in removedXscriptLst):
                                        lftvrXscriptLst.remove(xscript2)
                                        removedXscriptLst.append(xscript2)
                else:
                        continue
        
        # add leftovers as their own er set
        for leftover in lftvrXscriptLst:
                exRegSetLst.append([leftover])
        
        #DONE!!! all thats left is: configuring the data into proper table format in main
        #and ir stuff....
        return exRegSetLst

# change name
# description here
# "loops over all existing ER sets to add to existing set, otherwise creates new set in the form of a list"
def iterateAllERSet(xscript1, xscript2, erSetLst):
        for lst in erSetLst:
                for xscript in lst:
                        if xscript1 == xscript or xscript2 == xscript:
                                lst.append(xscript1)
                                lst.append(xscript2)
                                return;
                        else:
                                continue
                
        return [xscript1,xscript2]

def main():
        
        # input csv to dataframe
        inputDf = pd.read_csv(args.indir)
        
        return 'KSB'

if __name__ == '__main__':
        global args
        args = getOptions()
        
        inputDf = pd.read_csv(args.indir)
        
        if (args.includeIR == 'Y'):
                setTest = identifyERSet(inputDf, True)
        elif (args.includeIR == 'N'):
                setTest = identifyERSet(inputDf, False)
                
        print("complete")
        
        
        
        #output = pd.DataFrame.to_csv(outputDf)
        main()
        
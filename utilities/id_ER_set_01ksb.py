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
                print ("IR Retention excluded")
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
        
        # this could end up being useful?? it is useful!!
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
                                
                                #returns none if already in an existing set
                                if newLst is not None:
                                        exRegSetLst.append(newLst)
                                        
                        # removes transcript pair from leftover list (new variable name!!)
                        # removes if not already removed
                        # perhaps: move this to an already existing loop if possible, to save time
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
        return exRegSetLst, xscriptMstrLst

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

def createOutputDf(msterERSetLst, mstrXscriptLst):
        
        # ok so,,,, the output headers are the following: model_ID, set_num, num_pairs, num_in_set
        # model_ID = model_ID, extracted from each list in the master list
        # set_num = index + 1, index from master list
        # num_pairs = count of model_ID in that specific set
        # num_in_set = size of list, doesnt change between. ok. this shouldnt be bad
        
        # set up output dataframe w proper headers
        outDf = pd.DataFrame(columns=[
                                     'xscript_model_id', 
                                     'ER_set_num', 
                                     'xscript_freq_in_set', 
                                     'ER_set_size'])
        
        
        # add model IDs to dataframe
        outDf['xscript_model_id'] = mstrXscriptLst                
        
        # for every er set (list), we need to get:
                # the number of times a transcript appears
                # the number of transcripts in the lst
                #
        
        for er_set in msterERSetLst:
                break;
        
        return  outDf
        

def main():
        
        # input csv to dataframe
        inputDf = pd.read_csv(args.indir)
        # if (args.includeIR == 'Y'):
        #         setTest = identifyERSet(inputDf, True)
        #         outputTest = createOutputDf(setTest)
        # elif (args.includeIR == 'N'):
        #         setTest = identifyERSet(inputDf, False)
        #         outputTest = createOutputDf(setTest)

                
        # print("complete")
        
        return 'KSB'

# dont forget to move all of this to main()
if __name__ == '__main__':
        global args
        args = getOptions()
        
        inputDf = pd.read_csv(args.indir)
        
        if (args.includeIR == 'Y'):
                erSetLst, allXscriptLst = identifyERSet(inputDf, True)
                outputTest = createOutputDf(erSetLst, allXscriptLst)
        elif (args.includeIR == 'N'):
                erSetLst, allXscriptLst = identifyERSet(inputDf, False)
                outputTest = createOutputDf(erSetLst, allXscriptLst)

                
        print("complete")
        
        
        
        #output = pd.DataFrame.to_csv(outputDf)
        main()
        
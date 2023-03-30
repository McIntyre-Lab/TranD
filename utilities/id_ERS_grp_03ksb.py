#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:43:44 2023

@author: k.bankole
"""

"""
Identify possible exon region shared (ERS) groups using TRAND ouptput of a 1 or 2 GTF pairwise file 

Version 3 (Output 2 Files)

"""

import argparse
import pandas as pd
from collections import Counter
import os
import trand.io
import time

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        parser = argparse.ArgumentParser(description="Output a csv containing information on "
                                         "the exon region shared groups for a list of transcripts using TRAND output "
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
                choices=['Y','y','n','N'],
                help="Choose N to exclude transcript models with intron retention events from ERS groups.")
        
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



def idSharedExonRegion(inDf, intronRetention):
        """
        
        Take a dataframe from TranD output data, extract only the necesssary information
        and identify all of the ERS groups in the data. Outputs all of the ERS groups in list format as
        well as a list of all of the unique transcript IDs.

        Parameters
        ----------
        inDf : DATAFRAME
                TRAND OUTPUT DATA CONVERTED FROM CSV TO DATAFRAME.
        intronRetention : BOOLEAN
                DECIDES IF THE PRESENCE OF INTRON RETENTION WILL INCLUDE OR EXCLUDE XSCRIPTS FROM ERS GROUPS.

        Returns
        -------
        exRegShrdGrpLst : LIST (OF LISTS)
                LIST OF ALL ERS GROUPS (AS LISTS).
        unqXscriptLst : LIST (OF STRINGS)
                LIST OF ALL UNIQUE XSCRIPT IDs.

        """

        # Chop down input data frame into only the information on exon regions
        basicERInfoDf = inDf[
                [
                        "gene_id",
                        "transcript_1",
                        "transcript_2",
                        "prop_ER_similar",
                        "flag_IR"
                ]
        ].copy()
        
        # Concatenate Gene_ID onto transcript for later extraction
        basicERInfoDf['transcript_1'] = basicERInfoDf['gene_id'] + "/" + basicERInfoDf['transcript_1']
        basicERInfoDf['transcript_2'] = basicERInfoDf['gene_id'] + "/" + basicERInfoDf['transcript_2']

        # If excluding IR events, copy all IR flagged xscripts into a 
        # separate list for checking later (irSuspectLst)
        if (not intronRetention):                 

                irSspctLst = pd.concat([basicERInfoDf[basicERInfoDf['flag_IR']==1]['transcript_1'],
                                       basicERInfoDf[basicERInfoDf['flag_IR']==1]['transcript_2']]
                                       ).unique().copy().tolist()
        
                # Removes all IR involved xscripts from working dataframe
                erInfoDf = basicERInfoDf[basicERInfoDf["flag_IR"] == 0]
        
        # Otherwise do nothing to the data here, no xscripts involved in IR
        else: 
                print("IR Retention included")
                erInfoDf = basicERInfoDf
                irSspctLst = None



        # Chop IR info and gene_ID out of working Df (no longer necessary)
        erInfoDf = erInfoDf[
                [
                        'transcript_1',
                        'transcript_2',
                        'prop_ER_similar'
                ]
        ]        
                
        # DF REFORMATTING:
        # Convert ER similar proportion to boolean flag and...
        # "options" kills warning about chained assignment for pandas
        pd.options.mode.chained_assignment = None
        erInfoDf['prop_ER_similar'] = erInfoDf['prop_ER_similar'] == 1
        pd.options.mode.chained_assignment = "warn"
        
        #...rename the column accordingly (end reformatting)
        erInfoDf = erInfoDf.rename(columns={"prop_ER_similar":"flag_ER_full_overlap"})
        
        
        
        # Create the master list of all unique xscripts    
        unqXscriptLst = pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]).unique()
        
        # Create a list of leftover xscripts for checking that all xscripts are appear in an ERS group at least once
        lftvrLst = unqXscriptLst.copy().tolist()
        
        # Create a list to check against leftover xscripts to assure only correct xscripts are removed
        removedXscriptLst = []                
        
        # Convert Df to dictionary for faster iteration
        # Dictionary structure: 
                # key = column header
                # value = specific row value (xscript1, xscript2, er flag)
        # access flag_ER: row['flag_ER_full_overlap']
        # access xscript: row['transcript_1'] or row['transcript_2']
        iterDct = erInfoDf.to_dict('records')

        # Create empty master list of lists to contain all ERS groups        
        exRegShrdGrpLst = []  
        
        # Loop through every row in the inputDf
        for row in iterDct:
                fullOvlpFlag = row['flag_ER_full_overlap']
                xscript1 = row['transcript_1']
                xscript2 = row['transcript_2']
                
                # If there is an exon region with full overlap, do the following:
                if (fullOvlpFlag):
                        
                        # If its the first group (master list is empty), just create a new group with
                        # no looping 
                        if not exRegShrdGrpLst:
                                exRegShrdGrpLst.append([xscript1, xscript2])
                        
                        # Otherwise check all other lists to see if that xscript is already in an existing group
                        else:
                                newLst = checkAllERSGrps(xscript1=xscript1, 
                                                       xscript2=xscript2, 
                                                       ersGrpLst=exRegShrdGrpLst)
                                
                                # Above function call returns none if already in an existing group
                                if newLst is not None:
                                        exRegShrdGrpLst.append(newLst)
                                        
                        # Removes transcript pair from leftover list if not already removed
                        # Potential Improvement: move this to an already existing loop if possible, to save time
                        for lst in exRegShrdGrpLst:
                                if (xscript1 in lst) and (xscript1 not in removedXscriptLst):
                                        lftvrLst.remove(xscript1)
                                        if (not intronRetention and xscript1 in irSspctLst): irSspctLst.remove(xscript1)
                                        removedXscriptLst.append(xscript1)
                                if (xscript2 in lst) and (xscript2 not in removedXscriptLst):
                                        lftvrLst.remove(xscript2)
                                        if (not intronRetention and xscript2 in irSspctLst): irSspctLst.remove(xscript2)
                                        removedXscriptLst.append(xscript2)
                
                # Otherwise do nothing
                else:
                        continue
        
        # All scripts with no overlap should be in lftvrLst
        # Add leftovers as their own ERS group (they share with no one but themselves)
        for leftover in lftvrLst:
                exRegShrdGrpLst.append([leftover])
                if (not intronRetention and leftover in irSspctLst): irSspctLst.remove(leftover)

        
        # In above loops, have been removing any xscript added to a group from irSspct
        # Place all leftover irSspcts (every pair that the xscript is in has intron retention)
        # in their own group
        if irSspctLst:
                for irModel in irSspctLst:
                        exRegShrdGrpLst.append([irModel])        
        
        
        return exRegShrdGrpLst, unqXscriptLst

# i know there's an S but it really does not make sense if its not plural
def checkAllERSGrps(xscript1, xscript2, ersGrpLst):
        """
        Loops over all exisiting ERS groups to see if an xscript should be added to an existing group.
        Otherwise creates a new ERS group.

        Parameters
        ----------
        xscript1 : STRING
                ONE XSCRIPT ID TO CHECK.
        xscript2 : STRING
                OTHER XSCRIPT ID TO CHECK.
        ersGrpLst : LIST (OF LISTS)
                LIST OF ALL ALREADY EXISTING ERS GROUPS.

        Returns
        -------
        Quits loop and returns None if xscript already in existing group. 
        Otherwise returns a new LIST (new ERS group).
        
        """
        
        for lst in ersGrpLst:
                for xscript in lst:
                        
                        # Add xscript1 and 2 to group if already in an existing one, end loop.
                        if xscript1 == xscript or xscript2 == xscript:
                                lst.append(xscript1)
                                lst.append(xscript2)
                                return;
                                
                        # Otherwise keep looping
                        else:
                                continue
        
        # If looped through all groups and found nothing, create new group
        return [xscript1,xscript2]


def findNonolpPair(mstrERSGrpLst, anomalyXscript):
        """
        
        Helper function to create piped string containing a list of
        non-overlapping transcripts for the output dataframe

        Parameters
        ----------
        mstrERSGrpLst : LIST (OF LISTS)
                INPUT OF ALL ERS GROUPS IN LIST FORM (A LIST OF LISTS).
        anomalyXscript : STRING
                THE SPECIFIC XSCRIPT THAT HAS AT LEAST ONE NONOVERLAPPING PAIR
                WITHIN THE GROUP

        Returns
        -------
        outputString : STRING
                PIPED LIST OF TRANSCRIPTS IN STRING FORM TO BE OUTPUT TO DATAFRAME.

        """
        # Anomaly = Transcript that is in a group but does not have overlap 
        # with one of the xscripts in the group
        
        # Contains all unique Xscript Pairs that have the anomaly in them 
        # if the anomaly is on an odd index you want it and the one before
        # index = even or 0 you want it and the one after                        
        comparisonSet = set(())
        
        # Just a set of all unique transcripts in the group with the anomaly
        allXscriptSet = set(())       
        

        
        # look for the group that contains the anomaly
        for ers_grp in mstrERSGrpLst:
                
                # found it! 
                if anomalyXscript in ers_grp:
                        
                        # loop through the group containing the anomaly
                        for index, model in enumerate(ers_grp): 
                                
                                # grab just the xscript ID
                                xscript = model.split('/')[1]
                                
                                # creating list of all unique xscripts
                                allXscriptSet.add(xscript)
                                
                                # if at the anomaly, add it and the one
                                # before/after the comparison set
                                # essentially grabs the transcript pair
                                # makes sure to grab just the xscript (removes geneid)
                                if (model == anomalyXscript):
                                        comparisonSet.add(xscript)
                                        if index%2 == 0:
                                                comparisonSet.add(ers_grp[index+1].split('/')[1])
                                        else:
                                                comparisonSet.add(ers_grp[index-1].split('/')[1])
        
        # Create piped list of all xscripts that do not overlap with the anomaly
        outputString = "|".join(allXscriptSet - comparisonSet)
        
        return outputString


def createXscriptOutDf(mstrERSGrpLst, mstrXscriptLst):
        """
        Converts list output from idSharedExonRegion into a Dataframe
        
        Parameters
        ----------
        mstrERSGrpLst : LIST (OF LISTS)
                INPUT OF ALL ERS GROUPS IN LIST FORM (A LIST OF LISTS).
        mstrXscriptLst : LIST (OF STRINGS)
                INPUT OF ALL XSCRIPT IDs IN LIST FORM.

        Returns
        -------
        outDf : Dataframe
                OUTPUTS ALL ERS GROUPS IN DATAFRAME FORM USING DESCRIPTIVE HEADERS.

        """
           
        
        # Set up empty Df with proper headers
        outDf = pd.DataFrame(columns=[
                                     'gene_id',
                                     'xscript_model_id', 
                                     'ERS_grp_num', 
                                     'ERS_grp_size',
                                     'flag_nonolp_pair'])
        
        # Create a counter for each ERS Group for displaying xscript freq
        counterLst = [] # List of Counters for each ERS Group
        for ers_grp in mstrERSGrpLst:
                counterLst.append(dict(Counter(ers_grp)))
        
        # Note on the counter: it is a dictionary
        # key = xscript ID, value = the number of times the xscript appears in the group
        
        # Used for counting which group we are on, for displaying ERS Group num
        grpCount = 0;
                
        # Loop through all of the counters (1 for each group)
        for counter in counterLst:
                grpSize = len(counter)
                
                # Loop through all of the xscripts in each counter and add their info to the Df
                for key, value in counter.items():
                        
                        #Split xscript_id into an array, gene (0) and transcript (1)
                        keySplit = key.split('/')
                        geneID = keySplit[0]
                        xscriptID = keySplit[1]
                        
                        # flag if there is a pair with no overlap in the group
                        flag_nonolp = value < grpSize - 1
                        
                        # string containing possible list of nonoverlapping transcripts
                        nonolpXscript = None
                        
                        # if flag true, run function to create list of nonoverlapping transcripts
                        if flag_nonolp: 
                                nonolpXscript = findNonolpPair(mstrERSGrpLst=mstrERSGrpLst, anomalyXscript=key)
                        
                        # Create Row of DF With All Necessary Info
                        tmpDf = pd.DataFrame(
                                {'gene_id':[geneID],
                                 'xscript_model_id':[xscriptID], 
                                 'ERS_grp_num':[grpCount+1], 
                                 'ERS_grp_size':[grpSize],
                                 'flag_nonolp_pair':['1' if flag_nonolp else '0'],
                                 'nonolp_xscript_id':[nonolpXscript]
                                 })
                        
                        # Append Row to the Output Dataframe
                        outDf = pd.concat([outDf,tmpDf])

                grpCount+=1
        
        return outDf

def createERSOutDf(mstrERSGrpLst, mstrXscriptLst):
        
         
        # Set up empty Df with proper headers
        # Come up with ideas for other things to add to the columns
        outDf = pd.DataFrame(columns=[
                                     'gene_id',
                                     'ERS_grp_num', 
                                     'ERS_grp_size', 
                                     'flag_nonolp_pair',
                                     'flag_IR'])
        
        
        # move ERS_grp_size here
        # waga baga. boo.
        
        return outDf
        
        

# Run the Program
def main():
        
        # Input CSV to Df
        inputDf = pd.read_csv(args.indir)
        
        # Get input File Name
        input_file_name = os.path.splitext(os.path.basename(args.indir))[0]
        
        # Start timer to track how long the looping process takes
        tic = time.perf_counter()
        
        # Two options based on if IR is included or excluded
        if (args.includeIR.upper() == 'Y'):
                
                # List of all ERS Groups and all transcripts in the input data
                ersGrpLst, allXscriptLst = idSharedExonRegion(inDf=inputDf, intronRetention=True)
                
                # Converts above into 2 dfs to be output to csv
                xscriptDf = createXscriptOutDf(mstrERSGrpLst=ersGrpLst, mstrXscriptLst=allXscriptLst)
                ersDf = createERSOutDf(mstrERSGrpLst=ersGrpLst, mstrXscriptLst=allXscriptLst)

                # Configure descriptive file name
                xscript_output_file = "{}/{}_xscript_output.csv".format(args.outdir, input_file_name)
                ers_output_file = "{}/{}_ers_output.csv".format(args.outdir, input_file_name)
                
        elif (args.includeIR.upper() == 'N'):
                ersGrpLst, allXscriptLst = idSharedExonRegion(inDf=inputDf, intronRetention=False)
                
                xscriptDf = createXscriptOutDf(mstrERSGrpLst=ersGrpLst, mstrXscriptLst=allXscriptLst)
                ersDf = createERSOutDf(mstrERSGrpLst=ersGrpLst, mstrXscriptLst=allXscriptLst)
                
                xscript_output_file = "{}/{}_xscript_output_noIR.csv".format(args.outdir, input_file_name)
                ers_output_file = "{}/{}_ers_output.csv".format(args.outdir, input_file_name)


        # End timer to track how long the looping process takes
        toc = time.perf_counter()       
        print(f"complete, operation took {toc-tic:0.4f} seconds")
        
        # Assures that the outdir is empty or -f is enabled to overwrite existing directory
        trand.io.prepare_outdir(args.outdir, args.force)
        
        # Output Df to CSV
        xscriptDf.to_csv(xscript_output_file,index=False)
                        
        return ersGrpLst, xscriptDf

if __name__ == '__main__':
        global args
        args = getOptions()
        test, test2 = main()
        
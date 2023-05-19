#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:43:44 2023

@author: k.bankole
"""

"""
Identify possible exon region shared (ERS) groups using TRAND ouptput of a 1 or 2 GTF pairwise file 

Version 6: Fixes some bugs in the grouping logic (works...)
"""

import argparse
import pandas as pd
import numpy as np
import os
import time
import statistics as stats


class ERS_GRP:
        """
        Class representation of an ERS Group.
        
        num (int): The group number (main group identifier).
                
        gene_id (string): The gene that all the transcripts in the group belong to.
                
        num_er (int): Number of exon regions that all the transcripts in the group have.
                
        size (int): Number of transcripts in the group.
                
        xscriptSet (set of strings): A set of all the transcripts in the group (in string form).
                
        """
        
        def __init__(self, num, gene_id, num_er):
                self.num = num
                self.gene_id = gene_id
                self.num_er = num_er

                self.size = 0
                
                self.xscriptSet = set()
        
        # Add a transcript to the group
        def addXscript(self, xscript):
                self.xscriptSet.add(xscript)
                self.size+=1
        
        
        def __str__(self):
                return "ERS GROUP NUMBER: " + str(self.num) + ", XSCRIPT LIST: " + str(self.xscriptSet)
        
        def __iter__(self):
                return self
        
        def __next__(self):
                return self
        
        # Group numbers used to compare two groups.
        def __eq__(self, other):
                return self.num == other.num
        
        def __hash__(self):
                return hash(str(self))


class XSCRIPT:
        """
        Class representation of a transcript.
        
        xscript_id (string): The actual name of this transcript (main identifier).
                
        gene_id (string): The gene that this transcript belongs to.
                
        num_er (int): Number of exon regions that the transcript has.
                
        ovlpSet (set of strings): Set of all the transcripts that fully overlap with this transcript.
                
        ovlpCnt (int): Number of transcripts that this transcript has full overlap with.
                
        ers_grp_num (int): ERS group that the transcript belongs to.
                
        irSet (set of strings): Set of all the transcripts that overlap with this transcript and have intron retention activity.
                
        num_nuc_diff (list of ints): A list of all num NT diffs for all transcripts that have full overlap.
                
        prop_nuc_diff (list of floats): A list of all prop NT diffs for all transcripts that have full overlap.
        
        gtfOne (boolean): If the transcript from the first GTF (2 GTF input)
        
        gtfTwo (boolean): If the transcript from the second GTF (2 GTF input)
                
        """
        
        def __init__(self, xscript_id, gene_id, num_er):
                self.xscript_id = xscript_id
                self.gene_id = gene_id
                self.num_er = num_er

                self.ovlpSet = set()
                self.ovlpCnt = 0
                
                self.ers_grp_num = None
                
                self.irSet = set()
                
                self.num_nuc_diff = []
                self.prop_nuc_diff = []
                
                self.gtfOne = False
                self.gtfTwo = False
                
        # Add a transcript to the overlap set
        def addOlp(self, olp):
                self.ovlpSet.add(olp)
                self.ovlpCnt+=1
        
        # Add a transcript to the IR set
        def addIR(self, ir):
                self.irSet.add(ir)
        
        # Return whether or not this transcript has any IR activity
        def flagIR(self):
                return len(self.irSet) > 0
        
        # Add a num and prop diff to the list of nucleotide differences
        def addDiff(self, num, prop):
                self.num_nuc_diff.append(num)
                self.prop_nuc_diff.append(prop)
                
        def __eq__(self, other):
                return self.xscript_id == other.xscript_id
        
        def __str__(self):
                return "XSCRIPT: " + self.xscript_id + ", OVLPSET: " + str(self.ovlpSet) + ", NUM EXON REGIONS: " + str(self.num_er)
        
        
class GENE:
        """
        Class representation of a gene.
        
        gene_id (string): The name of the gene (main identifier).
                
        ersGrpSet (set of ERS_GRP objects): A set of all ERS groups belonging to this gene.
                
        numSet (int): The size of ersGrpSet.
                
        """
        
        def __init__(self, gene_id):
                self.gene_id = gene_id
                
                self.ersGrpSet = set()
                self.numGrp = 0
        
        def __str__(self):
                return self.gene_id
        
        def addGrp(self, grp):
                self.ersGrpSet.add(grp)
                self.numGrp += 1
        

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        parser = argparse.ArgumentParser(description="Output 3 csvs containing information on "
                                         "the exon region shared groups (one transcript focused, one gene focused, and one ERS group focused)"
                                         "for a list of transcripts using TRAND output "
                                         "data. Contains the option to include or exclude "
                                         "transcripts with intron retention events (--includeIR). "
                                         "Input a TRAND output file (csv) (--infile), "
                                         "IR inclusion option (--includeIR Y or N), and an output path (--outdir)."
                                         "Output directory must already exist. Also includes an option to"
                                         "use a prefix (--prefix) other than the original file name "
                                         "for the output files. "
                                         "Warning for 2 GTF files: Only genes with at least one "
                                         "transcript in both GTF files will be sorted into groups. "
                                         )
        
        # INPUT
        parser.add_argument(
                "-i",
                "--infile",
                dest="infile",
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
                help="Location of output directory, must already exist.")
        
        parser.add_argument(
                "-p",
                "--prefix",
                dest="prefix",
                help="Input a prefix for the output files. Defaults to original file name."
                )
        args = parser.parse_args()
        return args

def gleanInputDf(inDf, includeIR, gtfOne, gtfTwo):
        """
        
        Takes all the important information from the input dataframe. Converts
        all the information into a dictionary of XSCRIPT objects.

        Parameters
        ----------
        inDf : DATAFRAME
                Input dataframe from user input. CSV file that is converted to a dataframe.
                
        includeIR : BOOLEAN
                Whether or not transcripts with intron retention are included or excluded from ERS groups.

        Returns
        -------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.

        """
        
        # if ('transcript_in_gene' in inDf.columns):
                
        
        # Chop down input dataframe to only necessary information.
        erInfoDf = inDf[
                [
                        "gene_id",
                        "transcript_1",
                        "transcript_2",
                        "num_ER_T1_only",
                        "num_ER_T2_only",
                        "num_ER_shared",
                        "prop_ER_similar",
                        "num_nt_diff",
                        "prop_nt_diff",
                        "flag_IR"
                ]
        ].copy()
        
        #print (len(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]).unique()))

        # Convert all number rows into string for ease of access.
        erInfoDf['num_ER_shared'] = erInfoDf['num_ER_shared'].astype(str)
        erInfoDf['num_ER_T1_only'] = erInfoDf['num_ER_T1_only'].astype(str)
        erInfoDf['num_ER_T2_only'] = erInfoDf['num_ER_T2_only'].astype(str)
        erInfoDf['num_nt_diff'] = erInfoDf['num_nt_diff'].astype(str)
        erInfoDf['prop_nt_diff'] = erInfoDf['prop_nt_diff'].astype(str)
        
        unqXscriptSet = set(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]))
        print ("Number of transcripts: " + str(len(unqXscriptSet)))
        
        # Stick gene_id and number of exon region info onto transcript name (used for leftovers)
        if gtfOne and gtfTwo:
                erInfoDf['transcript_1'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_1'] + "/" + 
                        erInfoDf['num_ER_shared'].fillna(0) + "/" + 
                        erInfoDf['num_ER_T1_only'].fillna(0) + "/" +
                        gtfOne)
                
                erInfoDf['transcript_2'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_2'] + "/" + 
                        erInfoDf['num_ER_shared'].fillna(0) + "/" + 
                        erInfoDf['num_ER_T2_only'].fillna(0) + "/" +
                        gtfTwo)
        else:
                erInfoDf['transcript_1'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_1'] + "/" + 
                        erInfoDf['num_ER_shared'].fillna(0) + "/" + 
                        erInfoDf['num_ER_T1_only'].fillna(0))
                
                erInfoDf['transcript_2'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_2'] + "/" + 
                        erInfoDf['num_ER_shared'].fillna(0) + "/" + 
                        erInfoDf['num_ER_T2_only'].fillna(0)) 
        
        # Create set of all unique transcripts (used for leftovers)
        unqXscriptSet = set(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]))
        
        xscriptDct = {}
        
        # Create empty set used to check that transcripts have not already been added
        addedXscriptSet = set()
        
        # Convert working df to dictionary for faster iteration
        loopDct = erInfoDf.to_dict('records')
        
        # Loop through every row in the dataframe
        for row in loopDct:
                
                # Glean necessary info from the row and store it
                model1 = row['transcript_1']
                model2 = row['transcript_2']
                
                geneid = model1.split('/')[0]
                xscriptStr1 = model1.split('/')[1]
                xscriptStr2 = model2.split('/')[1]
                
                flagFullOvlp = row['prop_ER_similar'] == 1
                flagIR = row['flag_IR'] == 1
                
                numER = int(row['num_ER_shared'])
                
                numNTDiff = row['num_nt_diff']
                propNTDiff = row['prop_nt_diff']
                
                # Only perform further operations on the row if there is full overlap
                if (flagFullOvlp):
                        
                        # If the xscript has already been added, grab the already existing object from the dictionary
                        if (xscriptStr1 in addedXscriptSet):
                                tmpXscript1 = xscriptDct.get(xscriptStr1)
                                tmpXscript1.addDiff(numNTDiff, propNTDiff)
                        
                        # Otherwise, create a new XSCRIPT object 
                        else:
                                tmpXscript1 = XSCRIPT(xscript_id=xscriptStr1, gene_id=geneid, num_er=numER)
                                tmpXscript1.addDiff(numNTDiff, propNTDiff)
                        
                                if gtfOne:
                                        tmpXscript1.gtfOne = True
                        
                        # Also, the above adds the NT diff info to the XSCRIPT object once created
                        
                        
                        # Same for the other transcript in the pair        
                        if (xscriptStr2 in addedXscriptSet):
                                tmpXscript2 = xscriptDct.get(xscriptStr2)
                                tmpXscript2.addDiff(numNTDiff, propNTDiff)

                        else:
                                tmpXscript2 = XSCRIPT(xscript_id=xscriptStr2, gene_id=geneid, num_er=numER)
                                tmpXscript2.addDiff(numNTDiff, propNTDiff)
                                
                                if gtfTwo:
                                        tmpXscript2.gtfTwo = True
                                
                        
                        # If IR is included...
                        if (includeIR):
                                
                                # Full overlap = add each xscript to each others ovlpSet
                                tmpXscript1.addOlp(xscriptStr2)
                                tmpXscript2.addOlp(xscriptStr1)

                                # If there is IR as well, add each xscript to each others IR Set
                                if (flagIR):
                                        tmpXscript1.addIR(xscriptStr2)
                                        tmpXscript2.addIR(xscriptStr1)
                        
                        # If IR is excluded...
                        else:
                                
                                # Full overlap = add overlap only if there is no IR
                                # Again, if there is IR, add each xscript to each others IR Set
                                if (flagIR):
                                        tmpXscript1.addIR(xscriptStr2)
                                        tmpXscript2.addIR(xscriptStr1)
                                else:
                                        tmpXscript1.addOlp(xscriptStr2)
                                        tmpXscript2.addOlp(xscriptStr1)
                                        
                                        
                        # Add new objects to dictionary and "addedXscripts"
                        if (xscriptStr1 not in addedXscriptSet):
                                xscriptDct[xscriptStr1] = tmpXscript1
                                addedXscriptSet.add(xscriptStr1)

                        if (xscriptStr2 not in addedXscriptSet):
                                xscriptDct[xscriptStr2] = tmpXscript2
                                addedXscriptSet.add(xscriptStr2)
                
        
        # Creates a set of all transcripts that have full overlap
        olpXscriptSet = set(pd.concat([erInfoDf[erInfoDf['prop_ER_similar']==1]['transcript_1'],
                                 erInfoDf[erInfoDf['prop_ER_similar']==1]['transcript_2']]
                                ).unique())
        
        # Leftover transcripts that overlap with no other transcript:
        # All transcripts minus all transcripts that DO have overlap
        leftovers = unqXscriptSet - olpXscriptSet
        
        # Convert leftovers to XSCRIPT object and add to dictionary
        for leftover in leftovers:
                geneid = leftover.split('/')[0]
                xscriptStr = leftover.split('/')[1]
                
                if gtfOne and gtfTwo:
                        whichGTF = leftover.split('/')[4]
                
                
                if (xscriptStr not in addedXscriptSet):
                        # Number of exon regions = num_ER_shared + num_ER_T1_only
                        # or T2 only depending.
                        # This calculation works I promise. I think.
                        numER = int(leftover.split('/')[2]) + int(leftover.split('/')[3])
                        
                        tmpXscript = XSCRIPT(xscript_id=xscriptStr, gene_id=geneid, num_er=numER)
                        
                        if gtfOne and gtfTwo:
                                if whichGTF == gtfOne:
                                        tmpXscript.gtfOne = True
                                elif whichGTF == gtfTwo:
                                        tmpXscript.gtfTwo = True
                                else:
                                        print ("what happened")
                        
                        # NT diff does not matter for xscripts that will eventually be alone in a group
                        tmpXscript.addDiff(np.NaN,np.NaN)
                        
                        xscriptDct[xscriptStr] = tmpXscript
                        addedXscriptSet.add(xscriptStr)

        return xscriptDct, erInfoDf
        
# I know there's an S here but I really cannot think of a better name.
def createERSGrps(xscriptDct):                
        """
        
        Creates a list of ERS_GRP objects based on the XSCRIPT object dictionary.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.

        Returns
        -------
        ersGrpLst : LIST (OF ERS_GRPs)
                List of all ERS_GRPs.

        """
        
        # Create empty output list
        ersGrpLst = []
        #Create counter for ERS_group number
        groupCount = 0;
        
        # Loop through all XSCRIPT objects in the dictionary
        for xscript in xscriptDct.values():                                        
                # If there are any existing groups
                if ersGrpLst:
                        # Loop through each group
                        for ersGrp in ersGrpLst:
                                # If the transcript is in the set, add the transcript and all the transcript it overlaps with
                                # AND all the transcripts that overlap with what it overlaps with
                                if xscript.xscript_id in ersGrp.xscriptSet:
                                        ersGrp.addXscript(xscript.xscript_id)
                                        xscript.ers_grp_num = ersGrp.num                                                        
                                        ersGrp.xscriptSet.update(xscript.ovlpSet)
                                        
                                        ersGrp = addToGrp(grp=ersGrp, initXscript=xscript, xscriptDct=xscriptDct)
                                        
                                        break
                                
                        else:
                                # Otherwise just make a new group add the transcript (and all transcripts it overlaps with
                                # AND all the transcripts that overlap with what it overlaps with),
                                
                                groupCount += 1
                                tmpERS = ERS_GRP(num=groupCount, gene_id=xscript.gene_id, num_er=xscript.num_er)
                                
                                tmpERS.addXscript(xscript.xscript_id)
                                xscript.ers_grp_num = tmpERS.num
                                tmpERS.xscriptSet.update(xscript.ovlpSet)

                                tmpERS = addToGrp (grp=tmpERS, initXscript=xscript, xscriptDct=xscriptDct)
                                ersGrpLst.append(tmpERS) 
                                
                # If there are no groups yet, create group 1 and add the transcript (and all transcripts it overlaps with
                #AND all the transcripts that overlap with what it overlaps with),
                # add group to group list.
                else:
                        groupCount += 1
                        tmpERS = ERS_GRP(num=groupCount, gene_id=xscript.gene_id, num_er=xscript.num_er)
                        
                        tmpERS.addXscript(xscript.xscript_id)
                        xscript.ers_grp_num = tmpERS.num
                        tmpERS.xscriptSet.update(xscript.ovlpSet)
                        
                        tmpERS = addToGrp (grp=tmpERS, initXscript=xscript, xscriptDct=xscriptDct)
                        
                        ersGrpLst.append(tmpERS)                   
                        
        return ersGrpLst


# Adds all transcripts that overlap with what the transcripts overlaps with
def addToGrp(grp, initXscript, xscriptDct):
        """
        

        Parameters
        ----------
        grp : ERS_GRP
                Group to add an xscript to.
        initXscript : XSCRIPT
                Initial xscript to add to group
        xscriptDct : DICT (STR:XSCRIPT)
                The dictionary of all transcripts and xscript objects

        Returns
        -------
        grp : ERS_GRP
                Group with proper transcripts added.

        """
        
        grp.addXscript(initXscript.xscript_id)
        initXscript.ers_grp_num = grp.num
        grp.xscriptSet.update(initXscript.ovlpSet)
        
        
        # at this stage there may or may not be missing xscripts due to
        # the INITIAL xscript not overlapping with some other xscripts
        
        # here's the situation: 
                # transcript A overlaps with B, C, D
                # transcript B overlaps with A, D, E
                # initially, the group will contain A, B, C, D
                # and the group does not know about E, because A does not overlap with it
                # (but B does)
        
        # This is the solution to that (working on making it faster):
        smthAdded = True
        while smthAdded:
                
                smthAdded = False
                
                tmpSet = grp.xscriptSet
                
                for transcript in grp.xscriptSet:
                        tmpSize = len(tmpSet)
                        
                        tmpSet.update(xscriptDct.get(transcript).ovlpSet)
                        xscriptDct.get(transcript).ers_grp_num = grp.num
                        
                        if len(tmpSet) > tmpSize:
                                smthAdded = True
                                break
                
                grp.xscriptSet.update(tmpSet)
                grp.size = len(grp.xscriptSet)

        return grp
        
        
        
        
def createXscriptOutDf(xscriptDct, ersGrpLst):
        """
        
        Creates an output dataframe based on the information gleaned and stored
        in the XSCRIPT dictionary and ERS_GRP list. Transcript focused.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.
                
        ersGrpLst : LIST (OF ERS_GRPs)
                List of all ERS_GRPs.

        Returns
        -------
        outDf : DATAFRAME
                Output dataframe to be converted to a csv.

        """
        
        # Each column is given a list, each name represents a heading. -> Each position represents a row in the dataframe
        # Building these lists is far faster than building a dataframe for each row
        # and concatenating.
        geneIDLst = []
        xscriptStrLst = []
        grpNumLst = []
        nonOlpFlagLst = []
        nonOlpXscriptLst = []
        numERLst = []
        
        # Loop through every XSCRIPT object in the dictionary and append the necessary info to each list
        for xscript in xscriptDct.values():
                

                #Grab the ERS_GRP that the XSCRIPT belongs to for that information
                ersGrp = ersGrpLst[xscript.ers_grp_num - 1]
                
                grpSize = ersGrp.size
                
                # Determines whether there is a transcript that does not overlap with at least
                # one other transcript in the group
                flagNonOlp = xscript.ovlpCnt < grpSize - 1
                
                
                geneIDLst.append(xscript.gene_id)
                xscriptStrLst.append(xscript.xscript_id)
                grpNumLst.append(xscript.ers_grp_num)
                numERLst.append(xscript.num_er)
                
                nonOlpFlagLst.append('1' if flagNonOlp else '0')
                
                # If there is nonOlp. Create a piped list of all transcripts that the xscript does not overlap with
                if flagNonOlp:
                        grpXscriptSet = set()
                        grpXscriptSet.update(ersGrp.xscriptSet)
                        grpXscriptSet.remove(xscript.xscript_id)
                        
                        nonOlpXscript = "|".join(grpXscriptSet - xscript.ovlpSet)
                        nonOlpXscriptLst.append(nonOlpXscript)
                else:
                        nonOlpXscriptLst.append(np.NaN)

        # Create output dataframe using the lists
        outDf = pd.DataFrame(
                {
                'gene_id':geneIDLst,
                'xscript_model_id':xscriptStrLst, 
                'ERS_grp_num':grpNumLst, 
                'flag_nonolp_pair':nonOlpFlagLst,
                'nonolp_xscript_id':nonOlpXscriptLst,
                'num_ER':numERLst
                })
                
        return outDf

def createERSOutDf(ersGrpLst, xscriptDct, includeIR, gtfOne, gtfTwo):
        """
        Creates an output dataframe based on the information gleaned and stored
        in the XSCRIPT dictionary and ERS_GRP list. ERS_GRP focused.

        Parameters
        ----------
        ersGrpLst : LIST (OF ERS_GRPs)
                List of all ERS_GRPs.
                
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.
                
        includeIR : BOOLEAN
                Whether or not transcripts with intron retention are included or excluded from ERS groups.

        Returns
        -------
        outDf : DATAFRAME
                Output dataframe to be converted to a csv.

        """
        
        # Each column is given a list, each name represents a heading. -> Each position represents a row in the dataframe
        # Building these lists is far faster than building a dataframe for each row
        # and concatenating.
        numLst = []
        sizeLst = []
        geneIDLst = []
        xscriptStrLst = []
        nonOlpFlagLst = []
        irFlagLst = []
        numIRLst = []
        propIRLst = []
        numERLst = []
        minNumNTLst = []
        maxNumNTLst = []
        meanNumNTLst = []
        medNumNTLst = []
        minPropNTLst = []
        maxPropNTLst = []
        meanPropNTLst = []
        medPropNTLst = []
        sourceLst = []
        
        # Loop through every ERS_GRP in the list and append necessary info to each column list
        for ersGrp in ersGrpLst:
                numLst.append(ersGrp.num)
                sizeLst.append(ersGrp.size)
                geneIDLst.append(ersGrp.gene_id)
                numERLst.append(ersGrp.num_er)
                
                # Creat piped list of all the xscripts in the list
                xscriptStrLst.append("|".join(ersGrp.xscriptSet))
                
                # Create a list of all the num/prop NT diff of every transcript in the group
                numNTLst = []
                propNTLst = []
                containsGTFOne = False
                containsGTFTwo = False
                for xscriptStr in ersGrp.xscriptSet:        
                        
                        if (gtfOne and gtfTwo):
                                if xscriptDct.get(xscriptStr).gtfOne:
                                        containsGTFOne = True
                                
                                if xscriptDct.get(xscriptStr).gtfTwo:
                                        containsGTFTwo = True
                                

                        for numDiff in xscriptDct.get(xscriptStr).num_nuc_diff:
                                numNTLst.append(float(numDiff))
                                
                        for propDiff in xscriptDct.get(xscriptStr).prop_nuc_diff:
                                propNTLst.append(float(propDiff))
                
                if (gtfOne and gtfTwo):
                        if (containsGTFOne and containsGTFTwo):
                                sourceLst.append ("3")
                        elif (containsGTFOne):
                                sourceLst.append("1")
                        elif (containsGTFTwo):
                                sourceLst.append("2")
                        else:
                                print ("what happened")
                else:
                        sourceLst.append("0")
                
                # Do statistics on the list and add it to the appropriate columns
                minNumNTLst.append(min(numNTLst))
                maxNumNTLst.append(max(numNTLst))
                meanNumNTLst.append(stats.mean(numNTLst))
                medNumNTLst.append(stats.median(numNTLst))
                minPropNTLst.append(min(propNTLst))
                maxPropNTLst.append(max(propNTLst))
                meanPropNTLst.append(stats.mean(propNTLst))
                medPropNTLst.append(stats.median(propNTLst))
                
                # Find if there is at least one transcript with nonOlp
                for xscriptStr in ersGrp.xscriptSet:
                        if xscriptDct.get(xscriptStr).ovlpCnt < ersGrp.size - 1:
                                        nonOlpFlagLst.append('1')
                                        break
                else:
                                        nonOlpFlagLst.append('0')
                
                # Find the amount of xscripts with IR activity and set flag_ir accordingly
                if (includeIR):
                        irNumCnt = 0;
                        
                        for xscript in ersGrp.xscriptSet:
                                if xscriptDct.get(xscript).flagIR():
                                        irNumCnt += 1
                                        
                        if irNumCnt > 0:
                                irFlagLst.append('1')
                        else:
                                irFlagLst.append('0')
                                
                        numIRLst.append(irNumCnt)
                        propIRLst.append(irNumCnt/ersGrp.size)
                else:
                        irFlagLst.append('0')
                        numIRLst.append(0)
                        propIRLst.append(0/ersGrp.size)
        
        # Create output dataframe using the lists
        outDf = pd.DataFrame(
                {
                'ERS_grp_num':numLst, 
                'ERS_grp_size':sizeLst,
                'gene_id':geneIDLst,
                'xscripts':xscriptStrLst,
                'contains_which_gtf':sourceLst,
                'flag_nonolp_pair':nonOlpFlagLst,
                'flag_IR_in_set':irFlagLst,
                'num_IR_xscripts':numIRLst,
                'prop_IR':propIRLst,
                'num_ER':numERLst,
                'min_num_nt_diff':minNumNTLst,
                'max_num_nt_diff':maxNumNTLst,
                'mean_num_nt_diff':meanNumNTLst,
                'median_num_nt_diff':medNumNTLst,
                'min_prop_nt_diff':minPropNTLst,
                'max_prop_nt_diff':maxPropNTLst,
                'mean_prop_nt_diff':meanPropNTLst,
                'median_prop_nt_diff':medPropNTLst
                })
        
        return outDf

def createGeneOutDf(xscriptDct, ersGrpLst):
        """
        
        Creates a dictionary of GENE objects based on the information gleaned and stored
        in the XSCRIPT dictionary and ERS_GRP list. Then creates a gene focused 
        output dataframe.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.
                
        ersGrpLst : LIST (OF ERS_GRPs)
                List of all ERS_GRPs.

        Returns
        -------
        outDf : DATAFRAME
                Output dataframe to be converted to a csv.

        """
        
        # Create empty gene dictionary following the same structre as the XSCRIPT dictionary
        geneDct = {}
        
        # Loop through every group in the list.
        for grp in ersGrpLst:
                
                # If the gene_id that the group belongs to is already added,
                # Grab the already existing GENE object and add that group to its ersGrpSet
                if grp.gene_id in geneDct:
                        tmpGene = geneDct.get(grp.gene_id)
                        tmpGene.addGrp(grp)
                
                # Otherwise create a new GENE and add the group to its ersGrpSet, and add the 
                # new GENE to the dictionary
                else:
                        tmpGene = GENE(gene_id=grp.gene_id)
                        tmpGene.addGrp(grp)
                        
                        geneDct[grp.gene_id] = tmpGene
                        
        
        # Each column is given a list, each name represents a heading. -> Each position represents a row in the dataframe
        # Building these lists is far faster than building a dataframe for each row
        # and concatenating.
        geneIDLst = []
        numSetLst = []
        minERLst = []
        maxERLst = []
        meanERLst = []
        medERLst = []
        minSizeLst = []
        maxSizeLst = []
        meanSizeLst = []
        medSizeLst = []
        
        # Loop through every gene in the dictionary and append necessary info to each column list.
        for geneStr, gene in geneDct.items():
                
                geneIDLst.append(geneStr)
                numSetLst.append(len(gene.ersGrpSet))
                
                # Add every num_er and ers grp size for every group belonging to a gene to a list
                numERLst = []
                numSizeLst = []
                for grp in gene.ersGrpSet:
                        numERLst.append(int(grp.num_er))
                        numSizeLst.append(int(grp.size))
                
                if geneStr == "FBgn0039883":
                        print (numERLst)
                
                # and do stats on them.
                minERLst.append(min(numERLst))
                maxERLst.append(max(numERLst))
                meanERLst.append(stats.mean(numERLst))
                medERLst.append(stats.median(numERLst))
                
                minSizeLst.append(min(numSizeLst))
                maxSizeLst.append(max(numSizeLst))
                meanSizeLst.append(stats.mean(numSizeLst))
                medSizeLst.append(stats.median(numSizeLst))
        
        # Create output dataframe using the lists
        outDf = pd.DataFrame(
                {
                        'gene_id':geneIDLst,
                        'num_sets':numSetLst,
                        'min_ER':minERLst,
                        'max_ER':maxERLst,
                        'mean_ER':meanERLst,
                        'median_ER':medERLst,
                        'min_grp_size':minSizeLst,
                        'max_grp_size':maxSizeLst,
                        'mean_grp_size':meanSizeLst,
                        'median_grp_size':medSizeLst
                })
        
        return outDf  
                      
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


def main():
        """
        
        Runs the program.

        Raises
        ------
        OSError
                Raises an error if the output directory does not already exist.

        """
                
        
        # Get input File Name/Prefix
        if args.prefix:
                prefix = args.prefix
        else: 
                prefix = os.path.splitext(os.path.basename(args.infile))[0]
        
        # Used for testing, made it faster to run multiple times on the same file
        # if (os.path.exists(prefix + ".pickle") and os.path.getsize(prefix + ".pickle") > 0):
        #         with open (prefix + ".pickle", 'rb') as f:
        #                 inputDf = pickle.load(f)
        # else:
        #         inputDf = pd.read_csv (args.infile, low_memory=False)
        #         with open (prefix + ".pickle", 'wb') as f:
        #                 pickle.dump(inputDf, f)

        #Grab input DF from input CSV
        print("Reading input file.... ")
        inputDf = pd.read_csv (args.infile, low_memory=False)

        print ("Read complete! GTF Info: ")
        # Start timer to track how long the looping process takes
        tic = time.perf_counter()
        
        if (len(inputDf.columns) > 60):
                print("2 GTF")
                
                for column in inputDf.columns:
                        if column.startswith("num_transcript_in_gene_"):
                                gtfOne = column
                                break
                
                for column in inputDf.columns:
                        if column.startswith("num_transcript_in_gene_") and column != gtfOne:
                                gtfTwo = column
                                break 
                
                gtfOne = gtfOne[len(("num_transcript_in_gene_")):]
                gtfTwo = gtfTwo[len(("num_transcript_in_gene_")):]
                
                print ("GTF1: " + gtfOne)
                print ("GTF2: " + gtfTwo)
        else:
                print("1 GTF")
                gtfOne = None
                gtfTwo = None
        
        # Two options based on if IR is included or excluded
        if (args.includeIR.upper() == 'Y'):
                # Create XSCRIPT dictionary and ERS_GRP list
                mstrXscriptDct, erInfoDf = gleanInputDf(inDf=inputDf, includeIR=True, gtfOne=gtfOne, gtfTwo=gtfTwo)
                mstrERSGrpLst = createERSGrps(xscriptDct=mstrXscriptDct)
                
                
                # Create Output Dataframes
                xscriptDf = createXscriptOutDf(xscriptDct=mstrXscriptDct, ersGrpLst=mstrERSGrpLst)
                ersDf = createERSOutDf(ersGrpLst=mstrERSGrpLst, xscriptDct=mstrXscriptDct, includeIR=True, gtfOne=gtfOne, gtfTwo=gtfTwo)
                geneDf = createGeneOutDf(xscriptDct=mstrXscriptDct,ersGrpLst=mstrERSGrpLst)

                # Configure descriptive file name
                xscript_output_file = "{}/{}_xscript_output.csv".format(args.outdir, prefix)
                ers_output_file = "{}/{}_ers_output.csv".format(args.outdir, prefix)
                gene_output_file = "{}/{}_gene_output.csv".format(args.outdir, prefix)
                
                
        elif (args.includeIR.upper() == 'N'):  
                # Create XSCRIPT dictionary and ERS_GRP list
                mstrXscriptDct, erInfoDf = gleanInputDf(inDf=inputDf, includeIR=False, gtfOne=gtfOne, gtfTwo=gtfTwo)
                mstrERSGrpLst = createERSGrps(xscriptDct=mstrXscriptDct)
                
                # Create Output Dataframes
                xscriptDf = createXscriptOutDf(xscriptDct=mstrXscriptDct, ersGrpLst=mstrERSGrpLst)
                ersDf = createERSOutDf(ersGrpLst=mstrERSGrpLst, xscriptDct=mstrXscriptDct, includeIR=False, gtfOne=gtfOne, gtfTwo=gtfTwo)
                geneDf = createGeneOutDf(xscriptDct=mstrXscriptDct,ersGrpLst=mstrERSGrpLst)  

                # Configure descriptive file name
                xscript_output_file = "{}/{}_xscript_output_noIR.csv".format(args.outdir, prefix)
                ers_output_file = "{}/{}_ers_output_noIR.csv".format(args.outdir, prefix)
                gene_output_file = "{}/{}_gene_output_noIR.csv".format(args.outdir, prefix)

        
        # Output Df to CSV, will raise an error if output directory does not already exist
        try:
                xscriptDf.to_csv(xscript_output_file,index=False)
                ersDf.to_csv(ers_output_file,index=False)
                geneDf.to_csv(gene_output_file,index=False)
        except OSError:
                raise OSError("Output directory must already exist.")

        # End timer to track how long the process takes
        toc = time.perf_counter()       
        print(f"complete, operation took {toc-tic:0.4f} seconds")
        
        
        print("Generating split dataframe...")
        # Used in testing to see that each transcript is only assigned to one group.
        splitDf = split_column_by_sep(ersDf, col_name="xscripts", sep="|")
        dup = splitDf.duplicated("xscripts", keep=False)
        dupes = splitDf[dup]

        # splitDf.to_csv('test.csv')
        
        print ("Actually complete!")
        
        # Only returns these things for quick checking during development
        return mstrXscriptDct, mstrERSGrpLst, xscriptDf, ersDf, geneDf, erInfoDf, splitDf, dupes

if __name__ == '__main__':
        global args
        args = getOptions()
        mstrXscriptDct, mstrERSGrpLst, xscriptDf, ersDf, geneDf, erInfoDf, splitDf, dupes = main()

        

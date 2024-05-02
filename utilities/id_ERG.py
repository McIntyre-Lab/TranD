#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:43:44 2023

@author: k.bankole
"""

"""
Identify possible exon region groups (ERGs) using TRAND ouptput of a 1 or 2 GTF pairwise file 

Version 7/8: Possibly adding an output GTF file with info on each exon region for each group

Updated name from id_ERS_grp to id_ERG
"""

import argparse
import pandas as pd
import numpy as np
import os
import time
import statistics as stats
import trand.io
# import pickle


class ERG:
        """
        Class representation of an ERG.
        
        erg_id (int): The group identifier: gene_id + a number ranked by transcript (alphabetically).
                
        gene_id (string): The gene that all the transcripts in the group belong to.
                
        num_er (int): Number of exon regions that all the transcripts in the group have.
                
        size (int): Number of transcripts in the group.
                
        xscriptSet (set of strings): A set of all the transcripts in the group (in string form).
        
        exonChain (list of exons): Ordered list of exons
        """
        
        def __init__(self, num, gene_id, num_er):
                self.num = num
                self.gene_id = gene_id
                self.num_er = num_er

                self.size = 0
                
                self.erg_id = None
                
                self.xscriptSet = set()
                
                self.exonChain = []
                
                
        # Add a transcript to the group
        def addXscript(self, xscript):
                self.xscriptSet.add(xscript)
                self.size+=1
        
        def addExon(self, exon):
                self.exonChain.append(exon)
        
        def __str__(self):
                return "ERG NUMBER: " + str(self.num) + ", XSCRIPT LIST: " + str(self.xscriptSet)
        
        def __iter__(self):
                return self
        
        def __next__(self):
                return self
        
        # Group numbers used to compare two groups.
        def __eq__(self, other):
                return self.erg_id == other.erg_id
        
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
                
        erg_id (int): ERG that the transcript belongs to.
                
        irSet (set of strings): Set of all the transcripts that overlap with this transcript and have intron retention activity.
                
        num_nuc_noOvlp (list of ints): A list of all num NT diffs for all transcripts that have full overlap.
                
        prop_nuc_noOvlp (list of floats): A list of all prop NT diffs for all transcripts that have full overlap.
        
        gtfOne (boolean): If the transcript from the first GTF (2 GTF input)
        
        gtfTwo (boolean): If the transcript from the second GTF (2 GTF input)
        
        exonLst (list of EXONs): List of exons related* to the transcript
        
        *check documentation
        """
        
        def __init__(self, xscript_id, gene_id, num_er):
                self.xscript_id = xscript_id
                self.gene_id = gene_id
                self.num_er = num_er

                self.ovlpSet = set()
                self.ovlpCnt = 0
                
                self.erg_num = None
                self.erg_id = None
                
                self.irSet = set()
                
                self.num_nuc_noOvlp = []
                self.prop_nuc_noOvlp = []
                
                self.gtfOne = False
                self.gtfTwo = False
                
                self.exonLst = []
                
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
                self.num_nuc_noOvlp.append(num)
                self.prop_nuc_noOvlp.append(prop)
        
        def addExon(self, exon):
                if exon not in self.exonLst: 
                        self.exonLst.append(exon)
                
        def __eq__(self, other):
                return self.xscript_id == other.xscript_id
        
        def __str__(self):
                return "XSCRIPT: " + self.xscript_id + ", OVLPSET: " + str(self.ovlpSet) + ", NUM EXON REGIONS: " + str(self.num_er)
        
        
class GENE:
        """
        Class representation of a gene.
        
        gene_id (string): The name of the gene (main identifier).
                
        ergSet (set of ERG objects): A set of all ERGs belonging to this gene.
                
        numSet (int): The size of ergSet.
                
        """
        
        def __init__(self, gene_id):
                self.gene_id = gene_id
                
                self.ergSet = set()
                self.numGrp = 0
        
        def __str__(self):
                return self.gene_id
        
        def addGrp(self, grp):
                self.ergSet.add(grp)
                self.numGrp += 1
                
class EXON:
        """
        Class representation of an exon
        
        seqname (string): seqname
        
        start (int): start position of the exon 
        
        end (int): end position of the exon
        
        strand (int): strand
        """
        
        def __init__(self, seqname, start, end, strand):
                self.seqname = seqname
                self.start = start
                self.end = end
                self.strand = strand
                
        
        def __eq__(self, other):
                return self.start == other.start and self.end == other.end and self.seqname == other.seqname and self.strand == other.strand
        
        def __str__(self):
                return ("SEQNAME: " + self.seqname + ", START: " + str(self.start) + 
                        ", END: " + str(self.end) + ", STRAND: " + self.strand)
        

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        parser = argparse.ArgumentParser(description="Identifies groups of transcripts with full exon region overlap. "
                                         "Output 3 csvs containing information on "
                                         "these exon region groups (ERGs) (one transcript focused, one gene focused, and one ERG focused) "
                                         "for a list of transcripts using pairwise transcript distance TRAND output "
                                         "data. Contains the option to include or exclude "
                                         "transcripts with intron retention events (--includeIR). "
                                         "Input a pairwise transcript distance TRAND output file (csv) (--infile), "
                                         "IR inclusion option (--includeIR Y or N), and an output path (--outdir)."
                                         "Output directory must already exist. Also includes an option to"
                                         "use a prefix (--prefix) other than the original file name "
                                         "for the output files."
                                         "Allows the option to output a GTF file with representative"
                                         "exons for each group (--gtf)"
                                         "Warning for 2 GTF files: Only genes with at least one "
                                         "transcript in both GTF files will be sorted into groups. "
                                         )
        
        # INPUT
        parser.add_argument(
                "-i",
                "--infile",
                dest="infile",
                required=True,
                help="Location of pairwise transcript distance TRAND output file")
        
        parser.add_argument(
                "-ir",
                "--includeIR",
                dest="includeIR",
                required=True,
                default='Y',
                const='Y',
                nargs='?',
                choices=['Y','y','n','N'],
                help="Choose N to exclude transcript models with intron retention events from ERGs.")
        
        # OUTPUT
        parser.add_argument(
                "-o",
                "--outdir",
                dest="outdir",
                required=True,
                help="Location of output directory, must already exist.")
        
        parser.add_argument(
                "-g",
                "--gtf",
                dest="outputGTF",
                action="store_true",
                help="Output a GTF file with a representative transcript for each group."
                "Each exon region is represented by an exon in the GTF. Each exon is built"
                "using the earliest start and the latest end within the group."
                )
        
        parser.add_argument(
                "-w",
                "--which-gtf",
                dest="whichGTF",
                required=False,
                help="Add this command if entering a 1 GTF pairwise distance file. Sets the \"which_gtf\" "
                "column to whatever number is entered. Enter a 1 or a 2. Helpful if continuing on two compare this "
                "1 GTF ERG output to a 2 GTF output. Make sure that whichever number you enter matches "
                "the which_gtf in the other file. If 1 is melanogaster and 2 is simulans in the 2 GTF output, and you enter "
                "a 1 GTF file with melanogaster data, input a 1.")
        
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
                Whether or not transcripts with intron retention are included or excluded from ERGs.
                
        gtfOne: STRING
                Name of the first GTF (if 2GTF)
                =None if there is 1 GTF
                
        gtfTwo: BOOLEAN
                Name of the second GTF
                =None if there is 1 GTF


        Returns
        -------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.

        """
        
        # if ('transcript_inGene' in inDf.columns):

        #NEW
        # Chop down input dataframe to only necessary information.
        erInfoDf = inDf[
                [
                        "gene_id",
                        "transcript_1",
                        "transcript_2",
                        "num_ER_only_T1",
                        "num_ER_only_T2",
                        "num_ER_ovlp",
                        "prop_ER_ovlp",
                        "ER_only_T1",
                        "ER_only_T2",
                        "ER_ovlp",
                        "num_nt_noOvlp",
                        "prop_nt_noOvlp",
                        "flag_IR"
                ]
        ].copy()
        
        #print (len(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]).unique()))

        # Convert all number rows into string for ease of access.
        erInfoDf['num_ER_ovlp'] = erInfoDf['num_ER_ovlp'].astype(str)
        erInfoDf['num_ER_only_T1'] = erInfoDf['num_ER_only_T1'].astype(str)
        erInfoDf['num_ER_only_T2'] = erInfoDf['num_ER_only_T2'].astype(str)
        erInfoDf['num_nt_noOvlp'] = erInfoDf['num_nt_noOvlp'].astype(str)
        erInfoDf['prop_nt_noOvlp'] = erInfoDf['prop_nt_noOvlp'].astype(str)
        
        unqXscriptSet = set(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]))
        print ("Number of transcripts: " + str(len(unqXscriptSet)))
        print ("Number of genes: " + str(len(erInfoDf['gene_id'].unique())))
        
        
        # Stick gene_id and number of exon region info onto transcript name (used for leftovers)
        #NEW (add the junction string to each transcript)
        if gtfOne and gtfTwo:
                erInfoDf['transcript_1'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_1'] + "/" + 
                        erInfoDf['num_ER_ovlp'].fillna(0) + "/" + 
                        erInfoDf['num_ER_only_T1'].fillna(0) + "/" +
                        erInfoDf['ER_ovlp'].fillna('') + "/" +
                        erInfoDf['ER_only_T1'].fillna('') + "/" +
                        gtfOne)
                
                erInfoDf['transcript_2'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_2'] + "/" + 
                        erInfoDf['num_ER_ovlp'].fillna(0) + "/" + 
                        erInfoDf['num_ER_only_T2'].fillna(0) + "/" +
                        erInfoDf['ER_ovlp'].fillna('') + "/" +
                        erInfoDf['ER_only_T2'].fillna('') + "/" +
                        gtfTwo)
        else:
                erInfoDf['transcript_1'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_1'] + "/" + 
                        erInfoDf['num_ER_ovlp'].fillna(0) + "/" + 
                        erInfoDf['num_ER_only_T1'].fillna(0) + "/" +
                        erInfoDf['ER_ovlp'].fillna('') + "/" +
                        erInfoDf['ER_only_T1'].fillna(''))
                
                erInfoDf['transcript_2'] = (
                        erInfoDf['gene_id'] + "/" + 
                        erInfoDf['transcript_2'] + "/" + 
                        erInfoDf['num_ER_ovlp'].fillna(0) + "/" + 
                        erInfoDf['num_ER_only_T2'].fillna(0) + "/" +
                        erInfoDf['ER_ovlp'].fillna('') + "/" +
                        erInfoDf['ER_only_T2'].fillna(''))
                
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
                
                exonChain = model1.split('/')[4].split('|')
                
                flagFullOvlp = row['prop_ER_ovlp'] == 1
                flagIR = row['flag_IR'] == 1
                
                numER = int(row['num_ER_ovlp'])
                
                numNTDiff = row['num_nt_noOvlp']
                propNTDiff = row['prop_nt_noOvlp']
                
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
                                
                        
                        # exon stuff
                        for exon in exonChain:
                                seqname = exon.split(':')[0]
                                start = int(exon.split(':')[1])
                                end = int(exon.split(':')[2])
                                strand = exon.split(':')[3]
                                
                                tmpExon = EXON(seqname=seqname, start=start, end=end, strand=strand)
                                
                                tmpXscript1.addExon(tmpExon)
                                tmpXscript2.addExon(tmpExon)
                                        
                                        
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
        olpXscriptSet = set(pd.concat([erInfoDf[erInfoDf['prop_ER_ovlp']==1]['transcript_1'],
                                 erInfoDf[erInfoDf['prop_ER_ovlp']==1]['transcript_2']]
                                ).unique())
        
        # Leftover transcripts that overlap with no other transcript:
        # All transcripts minus all transcripts that DO have overlap
        leftovers = unqXscriptSet - olpXscriptSet
        
        # Convert leftovers to XSCRIPT object and add to dictionary
                # 0 = geneID
                # 1 = xscript
                # 2 = num_ER_ovlp
                # 3 = num_ER_only
                # 4 = ER_ovlp (junction string)
                # 5 = ER_only_T2 (junction string)
                # 6 = which GTF
        for leftover in leftovers:
                geneid = leftover.split('/')[0]
                xscriptStr = leftover.split('/')[1]
                
                sharedExonChain = leftover.split('/')[4].split('|')
                onlyExonChain = leftover.split('/')[5].split('|')
                
                exonChain = sharedExonChain + onlyExonChain
                                                
                if gtfOne and gtfTwo:
                        whichGTF = leftover.split('/')[6]
                
                if (xscriptStr not in addedXscriptSet):
                        # Number of exon regions = num_ER_ovlp + num_ER_only_T1
                        # or T2 only depending.
                        # This calculation works I promise. I think.
                        
                        numER = int(leftover.split('/')[2]) + int(leftover.split('/')[3])
                        
                        tmpXscript = XSCRIPT(xscript_id=xscriptStr, gene_id=geneid, num_er=numER)
                        
                        for exon in exonChain:
                                if exon != '':
                                        seqname = exon.split(':')[0]
                                        start = int(exon.split(':')[1])
                                        end = int(exon.split(':')[2])
                                        strand = exon.split(':')[3]
                                        
                                        tmpExon = EXON(seqname=seqname, start=start, end=end, strand=strand)
                                        tmpXscript.addExon(tmpExon)
                                
                        
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
def createERGs(xscriptDct):                
        """
        
        Creates a list of ERG objects based on the XSCRIPT object dictionary.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.

        Returns
        -------
        ergLst : LIST (OF ERGs)
                List of all ERGs.

        """
        
        # Create empty output list
        ergLst = []
        #Create counter for ERG number
        groupCount = 0;
        
        # Loop through all XSCRIPT objects in the dictionary
        for xscript in xscriptDct.values():
                
                # If there are any existing groups
                if ergLst:
                        # Loop through each group
                        for erg in ergLst:
                                # If the transcript is in the set, add the transcript and all the transcript it overlaps with
                                # AND all the transcripts that overlap with what it overlaps with
                                if xscript.xscript_id in erg.xscriptSet:
                                        erg.addXscript(xscript.xscript_id)

                                        xscript.erg_num = erg.num 

                                        erg.xscriptSet.update(xscript.ovlpSet)
                                        
                                        erg = addToGrp(grp=erg, initXscript=xscript, xscriptDct=xscriptDct)
                                        
                                        break
                                
                        else:
                                # Otherwise just make a new group add the transcript (and all transcripts it overlaps with
                                # AND all the transcripts that overlap with what it overlaps with),
                                
                                groupCount += 1
                                tmpERG = ERG(num=groupCount, gene_id=xscript.gene_id, num_er=xscript.num_er)
                                
                                tmpERG.addXscript(xscript.xscript_id)
                                
                                xscript.erg_num = tmpERG.num 
                                
                                tmpERG.xscriptSet.update(xscript.ovlpSet)

                                tmpERG = addToGrp (grp=tmpERG, initXscript=xscript, xscriptDct=xscriptDct)
                                ergLst.append(tmpERG) 
                                
                # If there are no groups yet, create group 1 and add the transcript (and all transcripts it overlaps with
                #AND all the transcripts that overlap with what it overlaps with),
                # add group to group list.
                else:
                        groupCount += 1
                        tmpERG = ERG(num=groupCount, gene_id=xscript.gene_id, num_er=xscript.num_er)
                        
                        xscript.erg_num = tmpERG.num 

                        tmpERG.xscriptSet.update(xscript.ovlpSet)
                        
                        tmpERG = addToGrp (grp=tmpERG, initXscript=xscript, xscriptDct=xscriptDct)
                        
                        ergLst.append(tmpERG)                   
                        
        return ergLst


# Adds all transcripts that overlap with what the transcripts overlaps with
def addToGrp(grp, initXscript, xscriptDct):
        """
        

        Parameters
        ----------
        grp : ERG
                Group to add an xscript to.
        initXscript : XSCRIPT
                Initial xscript to add to group
        xscriptDct : DICT (STR:XSCRIPT)
                The dictionary of all transcripts and xscript objects

        Returns
        -------
        grp : ERG
                Group with proper transcripts added.

        """
        
        # grp.addXscript(initXscript.xscript_id)
        # initXscript.erg_id = grp.erg_id
        # initXscript.erg_num = grp.num 

        # grp.xscriptSet.update(initXscript.ovlpSet)
        
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
                        xscriptDct.get(transcript).erg_num = grp.num
                        
                        if len(tmpSet) > tmpSize:
                                smthAdded = True
                                break
                
                grp.xscriptSet.update(tmpSet)
                grp.size = len(grp.xscriptSet)

        return grp

def createGeneDct(ergLst):
        # Create empty gene dictionary following the same structre as the XSCRIPT dictionary
        geneDct = {}
        
        # Loop through every group in the list.
        for grp in ergLst:
                
                # If the gene_id that the group belongs to is already added,
                # Grab the already existing GENE object and add that group to its ergSet
                if grp.gene_id in geneDct:
                        tmpGene = geneDct.get(grp.gene_id)
                        tmpGene.addGrp(grp)
                
                # Otherwise create a new GENE and add the group to its ergSet, and add the 
                # new GENE to the dictionary
                else:
                        tmpGene = GENE(gene_id=grp.gene_id)
                        tmpGene.addGrp(grp)
                        
                        geneDct[grp.gene_id] = tmpGene
        
        return geneDct

def idERGs(geneDct, xscriptDct):
        for geneStr, gene in geneDct.items():
                
                tmpERGLst = list(gene.ergSet)
                
                if len(tmpERGLst) > 1:
                        
                        for x in tmpERGLst:
                                x.xscriptSet = sorted(x.xscriptSet)
                                        
                        # tmpERGLst = sorted(tmpERGLst, key=lambda x: (len(x.xscriptSet), x))
                        
                        count = 1
                        for erg in tmpERGLst:
                                erg.erg_id = geneStr + "_" + str(count)
                                count += 1
                                
                                for xscript in erg.xscriptSet:
                                        xscriptDct.get(xscript).erg_id = erg.erg_id
                        
                        
                else:
                        tmpERGLst[0].erg_id = geneStr + "_1"
                        for xscript in tmpERGLst[0].xscriptSet:
                                xscriptDct.get(xscript).erg_id = tmpERGLst[0].erg_id
        return None
                
        
def createXscriptOutDf(xscriptDct, ergLst):
        """
        
        Creates an output dataframe based on the information gleaned and stored
        in the XSCRIPT dictionary and ERG list. Transcript focused.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.
                
        ergLst : LIST (OF ERGs)
                List of all ERGs.

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
        ergIDLst = []
        nonOlpFlagLst = []
        nonOlpXscriptLst = []
        numERLst = []
        whichGtfLst = []
        
        # Loop through every XSCRIPT object in the dictionary and append the necessary info to each list
        for xscript in xscriptDct.values():
                
                #Grab the ERG that the XSCRIPT belongs to for that information
                erg = ergLst[xscript.erg_num - 1]
                
                grpSize = erg.size
                
                # Determines whether there is a transcript that does not overlap with at least
                # one other transcript in the group
                flagNonOlp = xscript.ovlpCnt < grpSize - 1
                
                geneIDLst.append(xscript.gene_id)
                xscriptStrLst.append(xscript.xscript_id)
                ergIDLst.append(xscript.erg_id)
                numERLst.append(xscript.num_er)
                
                nonOlpFlagLst.append('1' if flagNonOlp else '0')
                
                
                if xscript.gtfOne and xscript.gtfTwo:
                        print("what")
                
                if xscript.gtfOne:
                        whichGtfLst.append(1)
                elif xscript.gtfTwo:
                        whichGtfLst.append(2)
                else:
                        whichGtfLst.append(None)
                
                # If there is nonOlp. Create a piped list of all transcripts that the xscript does not overlap with
                if flagNonOlp:
                        grpXscriptSet = set()
                        grpXscriptSet.update(erg.xscriptSet)
                        grpXscriptSet.remove(xscript.xscript_id)
                        
                        nonOlpXscript = "|".join(grpXscriptSet - xscript.ovlpSet)
                        nonOlpXscriptLst.append(nonOlpXscript)
                else:
                        nonOlpXscriptLst.append(np.NaN)
                        
        
        if all(x == None for x in whichGtfLst):
                whichGtfLst = args.whichGTF
        
                
        # Create output dataframe using the lists
        outDf = pd.DataFrame(
                {
                'gene_id':geneIDLst,
                'transcript_id':xscriptStrLst, 
                'ERG_id':ergIDLst, 
                'flag_nonOlp_pair':nonOlpFlagLst,
                'nonOlp_xscript_id':nonOlpXscriptLst,
                'num_ER':numERLst,
                })
        
        if whichGtfLst:
                outDf['which_gtf'] = whichGtfLst
        
                
        return outDf

def createERGOutDf(ergLst, xscriptDct, includeIR, gtfOne, gtfTwo):
        """
        Creates an output dataframe based on the information gleaned and stored
        in the XSCRIPT dictionary and ERG list. ERG focused.

        Parameters
        ----------
        ergLst : LIST (OF ERGs)
                List of all ERGs.
                
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.
                
        includeIR : BOOLEAN
                Whether or not transcripts with intron retention are included or excluded from ERGs.

        Returns
        -------
        outDf : DATAFRAME
                Output dataframe to be converted to a csv.

        """
        
        # Each column is given a list, each name represents a heading. -> Each position represents a row in the dataframe
        # Building these lists is far faster than building a dataframe for each row
        # and concatenating.
        ergIDLst = []
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
                
        
        # Loop through every ERG in the list and append necessary info to each column list
        for erg in ergLst:
                ergIDLst.append(erg.erg_id)
                sizeLst.append(erg.size)
                geneIDLst.append(erg.gene_id)
                numERLst.append(erg.num_er)
                
                # Creat piped list of all the xscripts in the list
                xscriptStrLst.append("|".join(erg.xscriptSet))
                
                # Create a list of all the num/prop NT diff of every transcript in the group
                numNTLst = []
                propNTLst = []
                containsGTFOne = False
                containsGTFTwo = False
                for xscriptStr in erg.xscriptSet:
                                       
                        if (gtfOne and gtfTwo):
                                if xscriptDct.get(xscriptStr).gtfOne:
                                        containsGTFOne = True
                                
                                if xscriptDct.get(xscriptStr).gtfTwo:
                                        containsGTFTwo = True
                                

                        for numDiff in xscriptDct.get(xscriptStr).num_nuc_noOvlp:
                                numNTLst.append(float(numDiff))
                                
                        for propDiff in xscriptDct.get(xscriptStr).prop_nuc_noOvlp:
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
                for xscriptStr in erg.xscriptSet:
                        if xscriptDct.get(xscriptStr).ovlpCnt < erg.size - 1:
                                        nonOlpFlagLst.append('1')
                                        break
                else:
                                        nonOlpFlagLst.append('0')
                
                # Find the amount of xscripts with IR activity and set flag_ir accordingly
                if (includeIR):
                        irNumCnt = 0;
                        
                        for xscript in erg.xscriptSet:
                                        
                                if xscriptDct.get(xscript).flagIR():
                                        irNumCnt += 1
                                        
                        if irNumCnt > 0:
                                irFlagLst.append('1')
                        else:
                                irFlagLst.append('0')
                                
                        numIRLst.append(irNumCnt)
                        propIRLst.append(irNumCnt/erg.size)
                else:
                        irFlagLst.append('0')
                        numIRLst.append(0)
                        propIRLst.append(0/erg.size)
        
        # Create output dataframe using the lists
        outDf = pd.DataFrame(
                {
                'ERG_id':ergIDLst, 
                'ERG_size':sizeLst,
                'gene_id':geneIDLst,
                'xscripts':xscriptStrLst,
                'contains_which_gtf':sourceLst,
                'flag_nonOlp_pair':nonOlpFlagLst,
                'flag_IR_in_set':irFlagLst,
                'num_IR_xscripts':numIRLst,
                'prop_IR':propIRLst,
                'num_ER':numERLst,
                'min_num_nt_noOvlp':minNumNTLst,
                'max_num_nt_noOvlp':maxNumNTLst,
                'mean_num_nt_noOvlp':meanNumNTLst,
                'median_num_nt_noOvlp':medNumNTLst,
                'min_prop_nt_noOvlp':minPropNTLst,
                'max_prop_nt_noOvlp':maxPropNTLst,
                'mean_prop_nt_noOvlp':meanPropNTLst,
                'median_prop_nt_noOvlp':medPropNTLst
                })
        
        return outDf

def createGeneOutDf(xscriptDct, ergLst, geneDct):
        """
        Creates a dictionary of GENE objects based on the information gleaned and stored
        in the XSCRIPT dictionary and ERG list. Then creates a gene focused 
        output dataframe.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the actual XSCRIPT object.
                
        ergLst : LIST (OF ERGs)
                List of all ERGs.

        Returns
        -------
        outDf : DATAFRAME
                Output dataframe to be converted to a csv.

        """
                        
        
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
                numSetLst.append(len(gene.ergSet))
                
                # Add every num_er and erg grp size for every group belonging to a gene to a list
                numERLst = []
                numSizeLst = []
                for grp in gene.ergSet:
                        numERLst.append(int(grp.num_er))
                        numSizeLst.append(int(grp.size))
                
                
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
                        'num_grps':numSetLst,
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

def buildGrpExonChain(ergLst, xscriptDct):
        """
        Builds the chain of exons for the representative "transcript" of each group.
        
        Parameters
        ----------
        ergLst : LIST (OF ERGs)
                List of all ERGs.
                
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the 
                actual XSCRIPT object.
            
        Returns
        -------
        Nothing.

        """
        
        for erg in ergLst:
                
                tmpExonLst = []
                
                for xscriptStr in erg.xscriptSet:
                        xscript = xscriptDct.get(xscriptStr)
                        
                        for exon in xscript.exonLst:
                                if exon not in tmpExonLst:
                                        tmpExonLst.append(exon)
                
                tmpExonLst.sort(key = lambda x: x.start)
                
                
                firstExon = True
                for thisExon in tmpExonLst:
                        if firstExon:
                                firstExon = False
                                
                                erg.addExon(thisExon)
                        else:
                                lastExon = erg.exonChain[-1]

                                if thisExon.start > lastExon.end:
                                        erg.addExon(thisExon)
                                else:
                                        if thisExon.end > lastExon.end:
                                                lastExon.end = thisExon.end
                                        else:
                                                continue  
                                
def createGTFDf(ergLst,xscriptDct):
        """
        Creates a dataframe for the exons within the representative transcript
        of each ERG. Compatible with trand.io.write_gtf.
        output dataframe.

        Parameters
        ----------
        xscriptDct : DICTIONARY (KEY = STR, VALUE = XSCRIPT)
                Dictionary of transcripts, with the key being the string and the value being the 
                actual XSCRIPT object.
                
        ergLst : LIST (OF ERGs)
                List of all ERGs.

        Returns
        -------
        exonDf : DATAFRAME
                Output dataframe to be converted to a GTF file ia trand.io.write_gtf()

        """
        
        buildGrpExonChain(ergLst=ergLst, xscriptDct=xscriptDct)
        
        seqnameLst = []
        startLst = []
        endLst = []
        strandLst = []
        ergIDLst = []
        geneIDLst = []

        for erg in ergLst:
                seqname = erg.exonChain[0].seqname
                strand = erg.exonChain[0].strand
                ergID = str(erg.erg_id)
                geneID = erg.gene_id
                
                
                for exon in erg.exonChain:
                        seqnameLst.append(seqname)
                        strandLst.append(strand)
                        ergIDLst.append(ergID)
                        geneIDLst.append(geneID)
                        
                        startLst.append(exon.start)
                        endLst.append(exon.end)
                
                        
        exonDf = pd.DataFrame(
                {
                        'seqname':seqnameLst,
                        'start':startLst,
                        'end':endLst,
                        'strand':strandLst,
                        'transcript_id':ergIDLst,
                        'gene_id':geneIDLst
                })
        
        return exonDf       

# Only used for dev testing
# def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
#         # Split variable by some character like '|' or ',' and keep all other values the same
#         if col_name == None:
#                 col_name = 'transcript_id'
                
#         if sep == None:
#                 sep = "|"
                
#         splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
#         splitList.index = splitList.index.droplevel(-1)
#         tempDF = df.copy()
#         del(tempDF[col_name])
#         splitDF = tempDF.join(splitList.rename(col_name))
#         if sort_list != None:
#                 splitDF = splitDF.sort_values(by=sort_list)
#         del(tempDF, splitList)
#         return splitDF


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
        
        #Used for testing, made it faster to run multiple times on the same file
       # if (os.path.exists(prefix + "csv.pickle") and os.path.getsize(prefix + "csv.pickle") > 0):
       #        with open (prefix + "csv.pickle", 'rb') as f:
       #                 inputDf = pickle.load(f)
       # else:
       #         inputDf = pd.read_csv (args.infile, low_memory=False)
       #         with open (prefix + "csv.pickle", 'wb') as f:
       #                 pickle.dump(inputDf, f)

        #Grab input DF from input CSV
        print("Reading input file.... ")
        
        inputDf = pd.read_csv (args.infile, low_memory=False)
        print ("Read complete! GTF Info: ")
        # Start timer to track how long the looping process takes
        omegatic = time.perf_counter()
        
        if (len(inputDf.columns) > 60):
                print("2 GTF")
                
                for column in inputDf.columns:
                        if column.startswith("num_transcript_inGene_"):
                                gtfOne = column
                                break
                
                for column in inputDf.columns:
                        if column.startswith("num_transcript_inGene_") and column != gtfOne:
                                gtfTwo = column
                                break 
                
                gtfOne = gtfOne[len(("num_transcript_inGene_")):]
                gtfTwo = gtfTwo[len(("num_transcript_inGene_")):]
                
                print ("GTF1: " + gtfOne)
                print ("GTF2: " + gtfTwo)
        else:
                print("1 GTF")
                gtfOne = None
                gtfTwo = None
        
        
        # Two options based on if IR is included or excluded
        
        if (args.includeIR.upper() == 'Y'):
                includeIR = True
        elif (args.includeIR.upper() == 'N'):
                includeIR = False
        else:
                print ("Input Error.", args.includeIR, "is not a valid IR option.")
                quit()

        
        # Create XSCRIPT dictionary and ERG list
        mstrXscriptDct, erInfoDf = gleanInputDf(inDf=inputDf, includeIR=includeIR, gtfOne=gtfOne, gtfTwo=gtfTwo)        
        
        toc = time.perf_counter()       
        print(f"Gleaning complete,  took {toc-omegatic:0.4f} seconds")
                
        tic = time.perf_counter()
        mstrERGLst = createERGs(xscriptDct=mstrXscriptDct)
        
        toc = time.perf_counter()       
        print(f"Grouping complete,  took {toc-tic:0.4f} seconds")
                
        
        tic = time.perf_counter()
        # Create Output Dataframes
        geneDct = createGeneDct(ergLst=mstrERGLst)
        idERGs(geneDct=geneDct, xscriptDct=mstrXscriptDct)
        
        ergDf = createERGOutDf(ergLst=mstrERGLst, xscriptDct=mstrXscriptDct, includeIR=includeIR, gtfOne=gtfOne, gtfTwo=gtfTwo)
        geneDf = createGeneOutDf(xscriptDct=mstrXscriptDct,ergLst=mstrERGLst, geneDct=geneDct)
        xscriptDf = createXscriptOutDf(xscriptDct=mstrXscriptDct, ergLst=mstrERGLst)

                
        # Configure descriptive file name
        xscript_output_file = "{}/{}_xscript_output.csv".format(args.outdir, prefix)
        erg_output_file = "{}/{}_erg_output.csv".format(args.outdir, prefix)
        gene_output_file = "{}/{}_gene_output.csv".format(args.outdir, prefix)
        
        gtf_output_file = "{}/{}_erg_gtf.gtf".format(args.outdir, prefix)
                               
        if not includeIR:  
                xscript_output_file = "{}/{}_xscript_output_noIR.csv".format(args.outdir, prefix)
                erg_output_file = "{}/{}_erg_output_noIR.csv".format(args.outdir, prefix)
                gene_output_file = "{}/{}_gene_output_noIR.csv".format(args.outdir, prefix)
                
                gtf_output_file = "{}/{}_erg_noIR_gtf.gtf".format(args.outdir, prefix)
        
        toc = time.perf_counter()       
        print(f"Df building complete,  took {toc-tic:0.4f} seconds")

        if args.outputGTF:
                gtfDf = createGTFDf(ergLst=mstrERGLst, xscriptDct=mstrXscriptDct)
        

        # Output Df to CSV, will raise an error if output directory does not already exist
        try:
                xscriptDf.to_csv(xscript_output_file,index=False)
                ergDf.to_csv(erg_output_file,index=False)
                geneDf.to_csv(gene_output_file,index=False)
                
                if args.outputGTF:
                        if os.path.isfile(gtf_output_file):
                                os.remove(gtf_output_file)
                        
                        trand.io.write_gtf(data=gtfDf, out_fhs={"gtf":gtf_output_file}, fh_name="gtf")
                
        except OSError:
                raise OSError("Output directory must already exist.")

        # End timer to track how long the process takes
        toc = time.perf_counter()       
        print(f"Complete! Operation took {toc-omegatic:0.4f} seconds")
        
        # Only returns these things for quick checking during development
        return mstrXscriptDct, mstrERGLst, xscriptDf, ergDf, geneDf, erInfoDf

if __name__ == '__main__':
        global args
        args = getOptions()
        mstrXscriptDct, mstrERGLst, xscriptDf, ergDf, geneDf, erInfoDf = main()

        

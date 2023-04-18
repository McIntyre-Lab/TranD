#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:43:44 2023

@author: k.bankole
"""

"""
Identify possible exon region shared (ERS) groups using TRAND ouptput of a 1 or 2 GTF pairwise file 

Version 5 (Turn ERS groups into an object!!)
"""

import argparse
import pandas as pd
import os
import trand.io
import time
from dataclasses import dataclass
from collections import Counter


# Using a class to make it far far far easier to add any further functionality
# to the utility if necessary

@dataclass
class ERS_GRP:
        # is num really necessary? they will be in a list
        # num_exon_regions: int
        
        def __init__(self, num, gene_id):
                self.num = num
                self.size = 0
                self.gene_id = gene_id
                self.xscriptSet = set()
        
        def addXscript(self, xscript):
                self.xscriptSet.add(xscript)
                self.size+=1
                
        def __str__(self):
                return "ERS GROUP NUMBER: " + str(self.num) + ", XSCRIPT LIST: " + str(self.xscriptSet)
        
        def __iter__(self):
                return self
        
        def __next__(self):
                return self


# explain each parameter
@dataclass
class XSCRIPT:        
        def __init__(self, xscript_id, gene_id):
                self.xscript_id=xscript_id
                self.gene_id=gene_id
                self.ovlpSet = set()
                self.ovlpCnt = 0
                self.ers_grp_num = None
                
                self.irSet = set()

                
        def __eq__(self, other):
                return self.xscript_id == other.xscript_id
        
        def __str__(self):
                return self.xscript_id + ", ovlp set: " + str(self.ovlpSet)
        
        def compare(self, other):
                return self.xscript_id == other
        
        def addOlp(self, olp):
                self.ovlpSet.add(olp)
                self.ovlpCnt+=1
                
        def addIR(self, ir):
                self.irSet.add(ir)
                
        def flagIR(self):
                return len(self.irSet) > 0

        
        # flag_IR: bool
        # num_exon_regions: int
        # num_nuc_diff: list
        # prop_nuc_diff: list

# i dont think this will be useful... lets keep it for now
@dataclass
class GENE:
        gene_id: str
        num_sets: int
        

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



def convertInputDataFrame(inDf, includeIR):
        
        erInfoDf = inDf[
                [
                        "gene_id",
                        "transcript_1",
                        "transcript_2",
                        "prop_ER_similar",
                        "flag_IR"
                ]
        ].copy()
        
        erInfoDf['transcript_1'] = erInfoDf['gene_id'] + "/" + erInfoDf['transcript_1']
        erInfoDf['transcript_2'] = erInfoDf['gene_id'] + "/" + erInfoDf['transcript_2']
        
        unqXscriptSet = set(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]))
        
        # basically. is there a way to search through the whole dataframe ONCE and pull out each xscript with all of its info
        # i think yes
        # i did it!!
                
        xscriptDct = {}
        addedXscripts = set()
        
        iterDct = erInfoDf.to_dict('records')
        for row in iterDct:
                # print()
                # print (row)
                # print()
                
                model1 = row['transcript_1']
                model2 = row['transcript_2']
                
                geneid = model1.split('/')[0]
                xscript1 = model1.split('/')[1]
                xscript2 = model2.split('/')[1]
                flagFullOvlp = row['prop_ER_similar'] == 1
                flagIR = row['flag_IR'] == 1
                
                
                if (flagFullOvlp):
                        
                        if (xscript1 in addedXscripts):
                                tmpXscript1 = xscriptDct.get(xscript1)
                        else:
                                tmpXscript1 = XSCRIPT(xscript1, geneid)
                                
                        if (xscript2 in addedXscripts):
                                tmpXscript2 = xscriptDct.get(xscript2)
                        else:
                                tmpXscript2 = XSCRIPT(xscript2, geneid)                                
                                
                                
                        # print ("xscript 2: " + xscript2)
                        # print ("xscript 1: " + xscript1)
                        # print ("oxscript 2: " + str(tmpXscript2))
                        # print ("oxscript 1: " + str(tmpXscript1))
                        
                        if (includeIR):
                                tmpXscript1.addOlp(xscript2)
                                tmpXscript2.addOlp(xscript1)
                                
                                if (flagIR):
                                        tmpXscript1.addIR(xscript2)
                                        tmpXscript2.addIR(xscript1)
                        else:
                                if (flagIR):
                                        tmpXscript1.addIR(xscript2)
                                        tmpXscript2.addIR(xscript1)
                                else:
                                        tmpXscript1.addOlp(xscript2)
                                        tmpXscript2.addOlp(xscript1)
                                        
                                        
                                                                
                
                        # print ("olp list 1: " + str(tmpXscript1.ovlpSet))
                        # print ("olp list 2: " + str(tmpXscript2.ovlpSet))
                        
                        
                        if (xscript1 not in addedXscripts):
                                xscriptDct[xscript1] = tmpXscript1
                                addedXscripts.add(xscript1)

                        if (xscript2 not in addedXscripts):
                                xscriptDct[xscript2] = tmpXscript2
                                addedXscripts.add(xscript2)
                
                
        allOlpXscripts = set(pd.concat([erInfoDf[erInfoDf['prop_ER_similar']==1]['transcript_1'],
                                 erInfoDf[erInfoDf['prop_ER_similar']==1]['transcript_2']]
                                ).unique())
        
        leftovers = unqXscriptSet - allOlpXscripts
        
        for leftover in leftovers:
                geneid = leftover.split('/')[0]
                xscript = leftover.split('/')[1]
                
                tmpXscript = XSCRIPT(xscript, geneid)
                xscriptDct[xscript] = tmpXscript
        
        # dont forget to deal with 2. xscripts completely removed due to IR(?)
        return xscriptDct

def xscriptToGrp(xscriptDct):                
        
        ersGrpLst = []
        
        groupCount = 0;
        for value in xscriptDct.values():
                # print ("a: ")
                # print (value)
                # print()
                
                if ersGrpLst:
                        # print ("oh rilly?")
                        # print ()
                        
                        for ersGrp in ersGrpLst:
                                
                                # print ("b: ")
                                # print (ersGrp)
                                # print()
                                
                                if value.xscript_id in ersGrp.xscriptSet:
                                        ersGrp.addXscript(value.xscript_id)
                                        value.ers_grp_num = ersGrp.num
                                        
                                        ersGrp.xscriptSet.update(value.ovlpSet)
                                        break
                        else:
                                groupCount += 1
                                tmpERS = ERS_GRP(groupCount, value.gene_id)
                                
                                tmpERS.addXscript(value.xscript_id)
                                value.ers_grp_num = tmpERS.num

                                tmpERS.xscriptSet.update(value.ovlpSet)
                                ersGrpLst.append(tmpERS)                                
                        
                else:
                        # print ("gene_id: ")
                        # print (value.gene_id)
                        groupCount += 1
                        tmpERS = ERS_GRP(groupCount, value.gene_id)
                        
                        # if i need values from the xscript (doubtful) i'll add __hash__ to xscript
                        tmpERS.addXscript(value.xscript_id)
                        value.ers_grp_num = tmpERS.num

                        tmpERS.xscriptSet.update(value.ovlpSet)

                        # print ("tmpERS: " + str(tmpERS))
                        # print()
                        ersGrpLst.append(tmpERS)                        
                        
        return ersGrpLst
        
        
def createXscriptOutDf(xscriptDct, ersGrpLst):
        # Set up empty Df with proper headers
        outDf = pd.DataFrame(columns=[
                                     'gene_id',
                                     'xscript_model_id', 
                                     'ERS_grp_num', 
                                     'flag_nonolp_pair'])
        
        geneIDLst = []
        xscriptIDLst = []
        grpNumLst = []
        nonOlpLst = []
        nonolpXscriptLst = []
        
        for xscript in xscriptDct.values():
                ersGrp = ersGrpLst[xscript.ers_grp_num - 1]
                
                grpSize = ersGrp.size
                flag_nonolp = xscript.ovlpCnt < grpSize - 1
                
                
                geneIDLst.append(xscript.gene_id)
                xscriptIDLst.append(xscript.xscript_id)
                grpNumLst.append(xscript.ers_grp_num)
                nonOlpLst.append('1' if flag_nonolp else '0')
                
                if flag_nonolp:
                        tmpSet = set()
                        tmpSet.update(ersGrp.xscriptSet)
                        tmpSet.remove(xscript.xscript_id)
                        
                        nonOlpXscript = "|".join(tmpSet - xscript.ovlpSet)
                        nonolpXscriptLst.append(nonOlpXscript)
                else:
                        nonolpXscriptLst.append(None)

        
        outDf = pd.DataFrame(
                {
                'gene_id':geneIDLst,
                'xscript_model_id':xscriptIDLst, 
                'ERS_grp_num':grpNumLst, 
                'flag_nonolp_pair':nonOlpLst,
                'nonolp_xscript_id':nonolpXscriptLst
                })
                
        return outDf

def createERSOutDf(ersGrpLst, xscriptDct, includeIR):
        outDf = pd.DataFrame(columns=[
                                     'ERS_grp_num', 
                                     'ERS_grp_size',
                                     'gene_id',
                                     'xscripts',
                                     'flag_nonolp_pair',
                                     'flag_IR_in_set',
                                     'num_IR_xscripts',
                                     'prop_IR'])
        
        numLst = []
        sizeLst = []
        geneIDLst = []
        xscriptLst = []
        flagNonOlpLst = []
        flagIRLst = []
        irNumLst = []
        propIRLst = []
        
        for ersGrp in ersGrpLst:
                numLst.append(ersGrp.num)
                sizeLst.append(ersGrp.size)
                geneIDLst.append(ersGrp.gene_id)
                xscriptLst.append("|".join(ersGrp.xscriptSet))
                
                for xscript in ersGrp.xscriptSet:
                        if xscriptDct.get(xscript).ovlpCnt < ersGrp.size - 1:
                                flagNonOlpLst.append('1')
                                break
                else:
                                flagNonOlpLst.append('0')
                
                
                if (includeIR):
                        for xscript in ersGrp.xscriptSet:
                                if xscriptDct.get(xscript).flagIR():
                                        flagIRLst.append('1')
                                        break
                        else:
                                        flagIRLst.append('0')
                        
                        irNumCnt = 0;
                        
                        for xscript in ersGrp.xscriptSet:
                                if xscriptDct.get(xscript).flagIR():
                                        irNumCnt += 1
                                
                        irNumLst.append(irNumCnt)
                        propIRLst.append(irNumCnt/ersGrp.size)
                else:
                        flagIRLst.append('0')
                        irNumLst.append(0)
                        propIRLst.append(0/ersGrp.size)

                
        outDf = pd.DataFrame(
                {
                'ERS_grp_num':numLst, 
                'ERS_grp_size':sizeLst,
                'gene_id':geneIDLst,
                'xscripts':xscriptLst,
                'flag_nonolp_pair':flagNonOlpLst,
                'flag_IR_in_set':flagIRLst,
                'num_IR_xscripts':irNumLst,
                'prop_IR':propIRLst
                })
        
        return outDf

def createGeneOutDf(xscriptDct, ersGrpLst): 
        outDf = pd.DataFrame(columns=[
                                     'gene_id', 
                                     'num_sets',
                                ])
        
        geneDct = Counter((grp.gene_id) for grp in ersGrpLst)
        
        outDf = pd.DataFrame(
                {
                        'gene_id':geneDct.keys(),
                        'num_sets':geneDct.values()
                })
        
        return outDf                        

        

def main():
        inputDf = pd.read_csv (args.indir)
        
        # Get input File Name
        input_file_name = os.path.splitext(os.path.basename(args.indir))[0]

        # Start timer to track how long the looping process takes
        tic = time.perf_counter()
        
        # Two options based on if IR is included or excluded
        if (args.includeIR.upper() == 'Y'):
                # Configure descriptive file name
                print("intron retention included")
                mstrXscriptDct = convertInputDataFrame(inDf=inputDf, includeIR=True)
                mstrERSGrpLst = xscriptToGrp(xscriptDct=mstrXscriptDct)

                xscript_output_file = "{}/{}_xscript_output.csv".format(args.outdir, input_file_name)
                ers_output_file = "{}/{}_ers_output.csv".format(args.outdir, input_file_name)
                gene_output_file = "{}/{}_gene_output.csv".format(args.outdir, input_file_name)
                
                ersDf = createERSOutDf(ersGrpLst=mstrERSGrpLst, xscriptDct=mstrXscriptDct, includeIR=True)
                
        elif (args.includeIR.upper() == 'N'):         
                mstrXscriptDct = convertInputDataFrame(inDf=inputDf, includeIR=False)
                mstrERSGrpLst = xscriptToGrp(xscriptDct=mstrXscriptDct)

                xscript_output_file = "{}/{}_xscript_output_noIR.csv".format(args.outdir, input_file_name)
                ers_output_file = "{}/{}_ers_output_noIR.csv".format(args.outdir, input_file_name)
                gene_output_file = "{}/{}_gene_output_noIR.csv".format(args.outdir, input_file_name)
                
                ersDf = createERSOutDf(ersGrpLst=mstrERSGrpLst, xscriptDct=mstrXscriptDct, includeIR=False)
                
        

        # Converts above into 2 dfs to be output to csv
        xscriptDf = createXscriptOutDf(xscriptDct=mstrXscriptDct, ersGrpLst=mstrERSGrpLst)
        geneDf = createGeneOutDf(xscriptDct=mstrXscriptDct,ersGrpLst=mstrERSGrpLst)


        # Output Df to CSV
        xscriptDf.to_csv(xscript_output_file,index=False)
        ersDf.to_csv(ers_output_file,index=False)
        geneDf.to_csv(gene_output_file,index=False)

        # End timer to track how long the looping process takes
        toc = time.perf_counter()       
        print(f"complete, operation took {toc-tic:0.4f} seconds")
        
        # Assures that the outdir is empty or -f is enabled to overwrite existing directory
        #trand.io.prepare_outdir(args.outdir, args.force)
                

        
                                
        return mstrXscriptDct, mstrERSGrpLst, xscriptDf, ersDf, geneDf

if __name__ == '__main__':
        global args
        args = getOptions()
        mstrXscriptDct, mstrERSGrpLst, xscriptDf, ersDf, geneDf = main()
        

def createGeneOutDf(mstrERSGrpLst, mstrXscriptLst, mstrGeneLst):
        
        outDf = pd.DataFrame(columns=[
                                     'gene_id', 
                                     'num_sets',
                                     'anotherheader'])


        
        for geneID in mstrGeneLst:
                
                setCount = 0;
                for ers_grp in mstrERSGrpLst:
                        
                        if (geneID == ers_grp[0].split("/")[0]):
                                setCount+=1
                        
                
                tmpDf = pd.DataFrame(
                        {'gene_id':[geneID],
                         'num_sets':[setCount],
                         'anotherheader':[0]
                        })
                outDf = pd.concat([outDf, tmpDf])
        return outDf


        
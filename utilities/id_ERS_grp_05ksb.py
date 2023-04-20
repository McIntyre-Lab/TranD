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
import numpy as np
import os
import time
from dataclasses import dataclass
import statistics as stats


# Using a class to make it far far far easier to add any further functionality
# to the utility if necessary

@dataclass
class ERS_GRP:
        
        def __init__(self, num, gene_id, num_er):
                self.num = num
                self.size = 0
                self.gene_id = gene_id
                self.xscriptSet = set()
                self.num_er = num_er
        
        def addXscript(self, xscript):
                self.xscriptSet.add(xscript)
                self.size+=1
                
        def __str__(self):
                return "ERS GROUP NUMBER: " + str(self.num) + ", XSCRIPT LIST: " + str(self.xscriptSet)
        
        def __iter__(self):
                return self
        
        def __next__(self):
                return self
        
        def __eq__(self, other):
                return self.num == other.num
        
        def __hash__(self):
                return hash(str(self))


# explain each parameter
@dataclass
class XSCRIPT:        
        def __init__(self, xscript_id, gene_id, num_er):
                self.xscript_id=xscript_id
                self.gene_id=gene_id
                self.ovlpSet = set()
                self.ovlpCnt = 0
                self.ers_grp_num = None
                self.num_er = num_er
                
                self.irSet = set()
                
                self.num_nuc_diff = []
                self.prop_nuc_diff = []

                
        def __eq__(self, other):
                return self.xscript_id == other.xscript_id
        
        def __str__(self):
                return self.xscript_id + ", ovlp set: " + str(self.ovlpSet) + ", exon regions: " + str(self.num_er)
        
        def addOlp(self, olp):
                self.ovlpSet.add(olp)
                self.ovlpCnt+=1
                
        def addIR(self, ir):
                self.irSet.add(ir)
                
        def flagIR(self):
                return len(self.irSet) > 0
        
        def addDiff(self, num, prop):
                self.num_nuc_diff.append(num)
                self.prop_nuc_diff.append(prop)

@dataclass
class GENE:
        
        def __init__(self, gene_id):
                self.gene_id = gene_id
                self.ersGrpSet = set()
                self.numSets = 0
        
        def __str__(self):
                return self.gene_id
        
        def addGrp(self, grp):
                self.ersGrpSet.add(grp)
                self.numSets += 1
        

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
                                         "IR inclusion option (--includeIR), and an output path (--outdir)."
                                         "Output directory must already exist. "
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
        
        args = parser.parse_args()
        return args

def convertInputDataFrame(inDf, includeIR):
        
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
        
        erInfoDf['num_ER_shared'] = erInfoDf['num_ER_shared'].astype(str)
        erInfoDf['num_ER_T1_only'] = erInfoDf['num_ER_T1_only'].astype(str)
        erInfoDf['num_ER_T2_only'] = erInfoDf['num_ER_T2_only'].astype(str)
        erInfoDf['num_nt_diff'] = erInfoDf['num_nt_diff'].astype(str)
        erInfoDf['prop_nt_diff'] = erInfoDf['prop_nt_diff'].astype(str)

        
        erInfoDf['transcript_1'] = (
                erInfoDf['gene_id'] + "/" + 
                erInfoDf['transcript_1'] + "/" + 
                erInfoDf['num_ER_shared'] + "/" + 
                erInfoDf['num_ER_T1_only'])
        
        erInfoDf['transcript_2'] = (
                erInfoDf['gene_id'] + "/" + 
                erInfoDf['transcript_2'] + "/" + 
                erInfoDf['num_ER_shared'] + "/" + 
                erInfoDf['num_ER_T2_only'])
        
        unqXscriptSet = set(pd.concat([erInfoDf['transcript_1'], erInfoDf['transcript_2']]))
                
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
                numER = row['num_ER_shared']
                numNTDiff = row['num_nt_diff']
                propNTDiff = row['prop_nt_diff']
                
                if (flagFullOvlp):
                        
                        if (xscript1 in addedXscripts):
                                tmpXscript1 = xscriptDct.get(xscript1)
                                tmpXscript1.addDiff(numNTDiff, propNTDiff)
                        else:
                                tmpXscript1 = XSCRIPT(xscript1, geneid, numER)
                                tmpXscript1.addDiff(numNTDiff, propNTDiff)

                                
                        if (xscript2 in addedXscripts):
                                tmpXscript2 = xscriptDct.get(xscript2)
                                tmpXscript2.addDiff(numNTDiff, propNTDiff)

                        else:
                                tmpXscript2 = XSCRIPT(xscript2, geneid, numER)
                                tmpXscript2.addDiff(numNTDiff, propNTDiff)
                                
                        
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
                  
                  
                  if (xscript not in addedXscripts): 
                          numER = int(leftover.split('/')[2]) + int(leftover.split('/')[3])
                          tmpXscript = XSCRIPT(xscript, geneid, numER)
                          tmpXscript.addDiff(np.NaN,np.NaN)
                          
                          xscriptDct[xscript] = tmpXscript
                          addedXscripts.add(xscript)
                          
        return xscriptDct

def xscriptToGrp(xscriptDct):                
        
        ersGrpLst = []
        
        groupCount = 0;
        for value in xscriptDct.values():
                
                if ersGrpLst:

                        for ersGrp in ersGrpLst:
                                
                                if value.xscript_id in ersGrp.xscriptSet:
                                        ersGrp.addXscript(value.xscript_id)
                                        value.ers_grp_num = ersGrp.num
                                        
                                        ersGrp.xscriptSet.update(value.ovlpSet)
                                        break
                        else:
                                groupCount += 1
                                tmpERS = ERS_GRP(groupCount, value.gene_id, value.num_er)
                                
                                tmpERS.addXscript(value.xscript_id)
                                value.ers_grp_num = tmpERS.num

                                tmpERS.xscriptSet.update(value.ovlpSet)
                                ersGrpLst.append(tmpERS)                                
                        
                else:
                        groupCount += 1
                        tmpERS = ERS_GRP(groupCount, value.gene_id, value.num_er)
                        
                        tmpERS.addXscript(value.xscript_id)
                        value.ers_grp_num = tmpERS.num

                        tmpERS.xscriptSet.update(value.ovlpSet)

                        ersGrpLst.append(tmpERS)                        
                        
        return ersGrpLst
        
        
def createXscriptOutDf(xscriptDct, ersGrpLst):        
        geneIDLst = []
        xscriptIDLst = []
        grpNumLst = []
        nonOlpLst = []
        nonolpXscriptLst = []
        numERLst = []
        
        for xscript in xscriptDct.values():
                ersGrp = ersGrpLst[xscript.ers_grp_num - 1]
                
                grpSize = ersGrp.size
                flag_nonolp = xscript.ovlpCnt < grpSize - 1
                
                
                geneIDLst.append(xscript.gene_id)
                xscriptIDLst.append(xscript.xscript_id)
                grpNumLst.append(xscript.ers_grp_num)
                numERLst.append(xscript.num_er)
                
                nonOlpLst.append('1' if flag_nonolp else '0')
                
                if flag_nonolp:
                        tmpSet = set()
                        tmpSet.update(ersGrp.xscriptSet)
                        tmpSet.remove(xscript.xscript_id)
                        
                        nonOlpXscript = "|".join(tmpSet - xscript.ovlpSet)
                        nonolpXscriptLst.append(nonOlpXscript)
                else:
                        nonolpXscriptLst.append(np.NaN)

        
        outDf = pd.DataFrame(
                {
                'gene_id':geneIDLst,
                'xscript_model_id':xscriptIDLst, 
                'ERS_grp_num':grpNumLst, 
                'flag_nonolp_pair':nonOlpLst,
                'nonolp_xscript_id':nonolpXscriptLst,
                'num_ER':numERLst
                })
                
        return outDf

def createERSOutDf(ersGrpLst, xscriptDct, includeIR):
        
        numLst = []
        sizeLst = []
        geneIDLst = []
        xscriptLst = []
        flagNonOlpLst = []
        flagIRLst = []
        irNumLst = []
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
        
        for ersGrp in ersGrpLst:
                numLst.append(ersGrp.num)
                sizeLst.append(ersGrp.size)
                geneIDLst.append(ersGrp.gene_id)
                numERLst.append(ersGrp.num_er)

                xscriptLst.append("|".join(ersGrp.xscriptSet))
                
                for xscript in ersGrp.xscriptSet:
                        
                        numNTLst = []
                        propNTLst = []
                        for numDiff in xscriptDct.get(xscript).num_nuc_diff:
                                numNTLst.append(float(numDiff))
                                
                        
                        for propDiff in xscriptDct.get(xscript).prop_nuc_diff:
                                propNTLst.append(float(propDiff))
                        
                        if xscriptDct.get(xscript).ovlpCnt < ersGrp.size - 1:
                                flagNonOlpLst.append('1')
                                break
                else:
                                flagNonOlpLst.append('0')
                
                
                minNumNTLst.append(min(numNTLst))
                maxNumNTLst.append(max(numNTLst))
                meanNumNTLst.append(stats.mean(numNTLst))
                medNumNTLst.append(stats.median(numNTLst))
                minPropNTLst.append(min(propNTLst))
                maxPropNTLst.append(max(propNTLst))
                meanPropNTLst.append(stats.mean(propNTLst))
                medPropNTLst.append(stats.median(propNTLst))
                
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
        geneDct = {}
        
        for grp in ersGrpLst:
                if grp.gene_id in geneDct:
                        tmpGene = geneDct.get(grp.gene_id)
                        tmpGene.addGrp(grp)
                else:
                        tmpGene = GENE(grp.gene_id)
                        tmpGene.addGrp(grp)
                        
                        geneDct[grp.gene_id] = tmpGene
                        
        
        geneIDLst = []
        numSetLst = []
        minERLst = []
        maxERLst = []
        meanERLst = []
        medERLst = []
        
        
        for key, value in geneDct.items():
                
                geneIDLst.append(key)
                numSetLst.append(len(value.ersGrpSet))
                
                numERLst = []
                for grp in value.ersGrpSet:
                        numERLst.append(int(grp.num_er))
                
                minERLst.append(min(numERLst))
                maxERLst.append(max(numERLst))
                meanERLst.append(stats.mean(numERLst))
                medERLst.append(stats.median(numERLst))
        
        outDf = pd.DataFrame(
                {
                        'gene_id':geneIDLst,
                        'num_sets':numSetLst,
                        'min_ER':minERLst,
                        'max_ER':maxERLst,
                        'mean_ER':meanERLst,
                        'median_ER':medERLst
                })
        
        return outDf                        

        

def main():
        inputDf = pd.read_csv (args.infile)
        
        # Get input File Name
        input_file_name = os.path.splitext(os.path.basename(args.infile))[0]

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
        try:
                xscriptDf.to_csv(xscript_output_file,index=False, encoding='utf-16')
                ersDf.to_csv(ers_output_file,index=False, encoding='utf-16')
                geneDf.to_csv(gene_output_file,index=False, encoding='utf-16')
        except OSError:
                raise OSError("Output directory must already exist. ")

        # End timer to track how long the looping process takes
        toc = time.perf_counter()       
        print(f"complete, operation took {toc-tic:0.4f} seconds")
        
        return mstrXscriptDct, mstrERSGrpLst, xscriptDf, ersDf, geneDf

if __name__ == '__main__':
        global args
        args = getOptions()
        mstrXscriptDct, mstrERSGrpLst, xscriptDf, ersDf, geneDf = main()


        
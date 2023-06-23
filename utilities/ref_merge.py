#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:50:48 2023

@author: k.bankole

Cross-Species Merge Key File

"""

import argparse
import pandas as pd
import sys

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(
            description=(
                    "Combine two maps created for different species and determine whether"
                    "transcripts between species \'match\' based on the following criteria: "
                    "1. Are RMP (which is a given due to how the map is created)"
                    "2. Are FSM OR ERG_noIR with a nt_noOvlp less than the number"
                    "input. Formula for determining nt_noOvlp given in the help for -s."
                    "If no -s is entered, transcripts are only matched if they are FSMs."
            )
        )
        
        parser.add_argument(
                "-1",
                "--map-1",
                dest="m1",
                required=True,
                help=(
                        "Input file location of the map created by "
                        "prelim_merge.py for the first reference."
                )
        )
        
        parser.add_argument(
                "-2",
                "--map-2",
                dest="m2",
                required=True,
                help=(
                        "Input file location of the map created by "
                        "prelim_merge.py for the second reference."
                )
        )
        
        parser.add_argument(
                "-n1",
                "--name1",
                dest="name1",
                required=True,
                help=(
                    "Name of the first reference (that matches prefix of the columns)"
                    "in the map."
                )
        )
        
        parser.add_argument(
                "-n2",
                "--name2",
                dest="name2",
                required=True,
                help=(
                    "Name of the second reference (that matches prefix of the columns)"
                    "in the map."
                )
        )
        
        parser.add_argument(
            "-s",
            "--small-nt-difference",
            dest="ntDiff",
            required=False,
            type=int,
            help="Threshold for the number of nucleotides to use to define a 'small' "
            "difference in ERG transcript pairs. " 
            "Calculate for any species by doing: (num nt per codon) * (average # exons per gene). Ex: "
            "Drosophila melanogaster, s = 3 * 5 = 15."
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

        

if __name__ == '__main__':
        # Parse command line arguments
        global args
        args = getOptions()
        
        inDf1 = pd.read_csv(args.m1, low_memory=False)
        inDf2 = pd.read_csv(args.m2, low_memory=False)
        
        print("Num Xscripts in Union 1: " + 
              str(len(set(pd.concat([inDf1['transcript_1'], inDf1['transcript_2']])))))

        print("Num Xscripts in Union 2: " + 
              str(len(set(pd.concat([inDf2['transcript_1'], inDf2['transcript_2']])))))
        
        if args.name1 + "_" not in inDf1.columns.any():
                sys.exit("Error: " + args.name1 + " is not found in any of the columns "
                         "in input 1. Please check that the prefix matches "
                         "the prefix on the columns in the input file.")
        
        if args.name2 + "_" not in inDf2.columns.any():
                sys.exit("Error: " + args.name2 + " is not found in any of the columns "
                         "in input 2. Please check that the prefix matches "
                         "the prefix on the columns in the input file.")
                
        
        union1 = pd.DataFrame()
        union2 = pd.DataFrame()

        union1['transcript_1'] = inDf1['transcript_1']
        union1['transcript_2'] = inDf1['transcript_2']
        
        union2['transcript_1'] = inDf2['transcript_1']
        union2['transcript_2'] = inDf2['transcript_2']
        
        #TODO: NRM???
                
        keepColLst =  ['gene_id', 
                       'num_nt_noOvlp',
                       'flag_FSM', 
                       'ERG_id',
                       'flag_ERG_match',
                       'flag_ERG_noIR',
                       'flag_ERG_wIR',
                       'flag_min_match',
                       'flag_RMP']

        for col in inDf1.columns:
                for keep in keepColLst:
                        if col[:len(args.name1 + "_" + keep)] == args.name1 + "_" + keep:
                                if "tie" not in col and "FSM_RMP" not in col:
                                        union1[col] = inDf1[col]
        
        
        for col in inDf2.columns:
                for keep in keepColLst:
                        if col[:len(args.name1 + "_" + keep)] == args.name2 + "_" + keep:
                                if "tie" not in col and "FSM_RMP" not in col:
                                        union2[col] = inDf2[col]        
        
        mergedDf = pd.merge(union1, union2, on=['transcript_1', 'transcript_2'], 
                         how='outer', indicator = 'merge_check')
        
        # print (mergedDf['merge_check'].value_counts(dropna=False).sort_index())

        oops = mergedDf[mergedDf['merge_check'].str.contains('right')]
        
        if len(oops) > 0:
                sys.exit ("Error, Union Mismatch. There are transcripts in one union that"
                          "do not appear in the other.")
        else:
                print ("Valid Merge. All transcripts appear in both inputs.")
        
        mergedDf.drop(columns="merge_check")
        
        print("Num Xscripts in Merged Df: " + 
              str(len(set(pd.concat([mergedDf['transcript_1'], mergedDf['transcript_2']])))))
        
        reordColLst = ['transcript_1', 'transcript_2']
        

        
        reorderedDf = pd.DataFrame()
        
        posCnt = 0
        for keep in keepColLst:
                if keep == 'ERG_id':
                        for col in mergedDf.columns:
                                if 'ERG_id' in col and col not in reordColLst:
                                        reordColLst.append(col)

                elif keep == 'flag_min_match':
                        for col in mergedDf.columns:
                                if 'flag_min_match' in col not in reordColLst:
                                        reordColLst.append(col)
                else:                        
                        reordColLst.append(args.name1 + "_" + keep)
                        reordColLst.append(args.name2 + "_" + keep)

        reorderedDf = mergedDf[reordColLst].copy()
        # TRANSCRIPT MATCH:
        # have to be RMP (duh)
        # is an fsm OR
        # is an ERS no IR w small nt diff
        ergMatch = ((reorderedDf[args.name1 + '_flag_ERG_noIR'] == 1) 
                            & (reorderedDf[args.name2 + "_flag_ERG_noIR"] == 1))
        
        fsmMatch = ((reorderedDf[args.name1 + "_flag_FSM"] == 1) 
                            & (reorderedDf[args.name2 + "_flag_FSM"] == 1))
        
        if args.ntDiff:
                ntDiffMatch = ((reorderedDf[args.name1 + "_num_nt_noOvlp"] < args.ntDiff) 
                                    & (reorderedDf[args.name2 + "_num_nt_noOvlp"] < args.ntDiff))
        
                transcriptMatch = fsmMatch | ergMatch & ntDiffMatch
        else:
                transcriptMatch = fsmMatch
        
        
        reorderedDf.insert(2, 'transcript_match', transcriptMatch * 1)
        
        outputFile = args.outDir + "/" + args.name1 + "_" + args.name2 + "_species_merge.csv"
        
        try:
                reorderedDf.to_csv(outputFile,index=False)
        except OSError:
                raise OSError("Output directory must already exist.")
        
        print ("Complete!")                
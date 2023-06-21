#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:50:48 2023

@author: k.bankole
"""

import argparse
import pandas as pd

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(
            description=(
                    "Combine two unions created on different references."
            )
        )
        
        parser.add_argument(
                "-u1",
                "--union-1",
                dest="u1",
                required=True,
                help=(
                        "Input file location of the union created by"
                        "prelim_merge.py for the first reference."
                )
        )
        
        parser.add_argument(
                "-u2",
                "--union-2",
                dest="u2",
                required=True,
                help=(
                        "Input file location of the union created by"
                        "prelim_merge.py for the second reference."
                )
        )
        
        parser.add_argument(
                "-n1",
                "--name1",
                dest="name1",
                required=True,
                help=(
                    "Name of the second reference (that matches prefix of the columns)"
                    "In the union."
                )
        )
        
        parser.add_argument(
                "-n2",
                "--name2",
                dest="name2",
                required=True,
                help=(
                    "Name of the second reference (that matches prefix of the columns)"
                    "In the union."
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

        

if __name__ == '__main__':
        # Parse command line arguments
        global args
        args = getOptions()
        
        inDf1 = pd.read_csv(args.u1, low_memory=False)
        inDf2 = pd.read_csv(args.u2, low_memory=False)
        
        
        print("Num Xscripts in Union 1: " + 
              str(len(set(pd.concat([inDf1['transcript_1'], inDf1['transcript_2']])))))

        print("Num Xscripts in Union 2: " + 
              str(len(set(pd.concat([inDf2['transcript_1'], inDf2['transcript_2']])))))
        
        union1 = inDf1[
                [
                        args.name1 + "_gene_id",
                        "transcript_1",
                        "transcript_2",
                        args.name1 + "_num_nt_noOvlp",
                        args.name1 + "_flag_FSM",
                        args.name1 + "_flag_ERG_match",
                ]
        ].copy() 
        
        for col in inDf1.columns:
                if "ERG_id" in col:
                        union1[col] = inDf1[col]
        
        union2 = inDf2[
                [
                        args.name2 + "_gene_id",
                        "transcript_1",
                        "transcript_2",
                        args.name2 + "_num_nt_noOvlp",
                        args.name2 + "_flag_FSM",
                        args.name2 + "_flag_ERG_match",
                ]
        ].copy() 
        
        for col in inDf2.columns:
                if "ERG_id" in col:
                        union2[col] = inDf2[col]
        
        

        
        
        mergedDf = pd.merge(union1, union2, on=['transcript_1', 'transcript_2'], 
                         how='outer', indicator = 'merge_check')
        
        # print (mergedDf['merge_check'].value_counts(dropna=False).sort_index())

        oops = mergedDf[mergedDf['merge_check'].str.contains('right')]
        
        if len(oops) > 0:
                print ("Error, Union Mismatch")
        else:
                print ("Valid Merge")
        
        mergedDf.drop(columns="merge_check")
        reorderedDf = mergedDf[[
                        "transcript_1",
                        "transcript_2",
                        args.name1 + "_gene_id", 
                        args.name2 + "_gene_id", 
                        args.name1 + "_num_nt_noOvlp",
                        args.name2 + "_num_nt_noOvlp",
                        args.name1 + "_flag_FSM",
                        args.name2 + "_flag_FSM",
                        args.name1 + "_flag_ERG_match",
                        args.name2 + "_flag_ERG_match",
                        
                ]
        ].copy()
        
        for col in mergedDf.columns:
                if "ERG_id" in col:
                        reorderedDf[col] = mergedDf[col]
        
        print ("Transcripts that are the \"same:\"")
        for row in mergedDf.to_dict('records'):
                if row[args.name1 + "_flag_ERG_match"] == row[args.name2 + "_flag_ERG_match"]:
                        if row[args.name1 + "_flag_FSM"] == row[args.name2 + "_flag_FSM"]:
                                if row[args.name1 + "_num_nt_noOvlp"] < 15 and row[args.name2 + "_num_nt_noOvlp"] < 15: 
                                        print (row['transcript_1'], row['transcript_2'])


        print ("Complete!")                
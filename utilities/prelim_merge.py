#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:04:12 2023

@author: k.bankole

This is the TMM (ERG + PD)
"""

import argparse
import pandas as pd

def getOptions():
        # Parse command line arguments
        parser = argparse.ArgumentParser(
            description=(
                "Create a Transcript Map File from TranD output files."
            )
        )
        
        parser.add_argument(
                "-e",
                "--ERG-file",
                dest="ERG",
                required=True,
                help=(
                    "Input file location for xscript_output ERG file created from TranD output data."
                )
        )
        
        parser.add_argument(
                "-p",
                "--pairwise",
                dest="PD",
                required=True,
                help=(
                    "Input file location of pairwise distance file."
                )
        )
        
        parser.add_argument(
                "-1",
                "--name1",
                dest="GTF1",
                required=False,
                default="gtf1",
                help=(
                    "Optional name for the first GTF when processing 2 GTF output."
                    "Added to the header of the column when merged."
                    "Default: gtf1"
                )
        )
        
        parser.add_argument(
                "-2",
                "--name2",
                dest="GTF2",
                required=False,
                default="gtf2",
                help=(
                    "Optional name for the first GTF when processing 2 GTF output."
                    "Added to the header of the column when merged."
                    "Default: gtf2"
                )
        )
                
        parser.add_argument(
                "-o",
                "--outdir",
                dest="outDir",
                required=True,
                help=(
                        "Output directory, must already exist."
                )
        )
        
        parser.add_argument(
                "-c",
                "--col-prefix",
                dest="colPrefix",
                required=False,
                default = None,
                help = (
                        "Add a short prefix for all the columns"
                        "in the output. Necessary if comparing species."
                        "Default: no prefix"
                )
        )
        
        args = parser.parse_args()
        return args

        

if __name__ == '__main__':
        # Parse command line arguments
        global args
        args = getOptions()
        
        ergDf = pd.read_csv(args.ERG, low_memory=False)
        print ("ERG INFO: ")
        print ("ERG XSCRIPTS: " + str(len(ergDf)))
        
        ergDf1 = ergDf[ergDf['which_gtf'] == 1]
        print ("GTF1 XSCRIPTS: " + str(len(ergDf1)))
        
        ergDf2 = ergDf[ergDf['which_gtf'] == 2]
        print ("GTF2 XSCRIPTS: " + str(len(ergDf2)))
        
        print ("GENES: " + str(len(set(ergDf['gene_id']))))
        
        ergDf1 = ergDf1.rename(columns = {"xscript_model_id":"transcript_1", 
                                 "ERG_id":"ERG_id_" + args.GTF1, 
                                 "flag_nonOlp_pair":"flag_nonOlp_pair_" + args.GTF1,
                                 "nonOlp_xscript_id":"nonOlp_xscript_id_" + args.GTF1,
                                 "num_ER":"num_ER_" + args.GTF1})
        
        ergDf2 = ergDf2.rename(columns = {"xscript_model_id":"transcript_2", 
                                 "ERG_id":"ERG_id_" + args.GTF2, 
                                 "flag_nonOlp_pair":"flag_nonOlp_pair_" + args.GTF2,
                                 "nonOlp_xscript_id":"nonOlp_xscript_id_" + args.GTF2,
                                 "num_ER":"num_ER_" + args.GTF2})
        
        
        print ()
        print ("PAIRWISE DISTANCE INFO: ")        

        pdDf = pd.read_csv(args.PD, low_memory=False)
        
        print ("Len PD: " + str(len(pdDf)))
        print ("Genes: " + str(len(set(pdDf['gene_id']))))
        print ("Len T1: " + str(len(set(pdDf['transcript_1']))))
        print ("Len T2: " + str(len(set(pdDf['transcript_2']))))

        unqXscriptSet = set(pd.concat([pdDf['transcript_1'], pdDf['transcript_2']]))
        
        print()
        print ("Xscripts: " + str(len(unqXscriptSet)))        
        
        print()
        
        
        for col in pdDf.columns:
                
                if "flag_min_match_" in col:
                        gtf1 = col[len(("flag_min_match_")):]
                        break
        
        for col in pdDf.columns:
                if "flag_min_match_" in col:
                        gtf2 = col[len(("flag_min_match_")):]
                        
                        if gtf2 == gtf1:
                                continue
                        else:
                                break
        
        
        minmatchDf = pdDf[(pdDf['flag_min_match_' + gtf1] == 1) | (pdDf['flag_min_match_' + gtf2] == 1)]
        print ()
        print ("flagminmatch: " + str(len(minmatchDf)))
        
        trueMinDf = minmatchDf[minmatchDf['flag_RMP'] == 1]
        print ("Num true minimums: " + str(len(trueMinDf)))
        print ()
        
        
        # Merging Time
        ergDf1 = ergDf1.drop(columns="gene_id")
        ergDf2 = ergDf2.drop(columns="gene_id")
        ergDf1 = ergDf1.drop(columns="which_gtf")
        ergDf2 = ergDf2.drop(columns="which_gtf")
        
        # First Merge
        merge1 = pd.merge(trueMinDf, ergDf1, on=['transcript_1'], how='outer', indicator='merge_check')
        
        print ("First Merge Check: ")
        print (merge1['merge_check'].value_counts(dropna=False).sort_index())
        print()
        
        oops1 = merge1[merge1['merge_check'].str.contains('right')]
        
        # Remove Xscripts that dont appear in the subset PD
        subMerge1 = merge1[~merge1['merge_check'].str.contains('right')]
        subMerge1 = subMerge1.drop(columns="merge_check")
        
        # Second Merge
        merge2 = pd.merge(subMerge1, ergDf2, on=['transcript_2'], how='outer', indicator = 'merge_check')
        
        print ("Second Merge Check: ")
        print (merge2['merge_check'].value_counts(dropna=False).sort_index())
        print()
        
        oops2 = merge2[merge2['merge_check'].str.contains('right')]
        
        # Remove Xscripts that dont appear in the subset PD
        # Union Complete
        unionDf = merge2[~merge2['merge_check'].str.contains('right')].copy()


        #Checking Removed Transcripts
        # Create List of Removed Transcipts
        minXscripts = set(pd.concat([trueMinDf['transcript_1'], trueMinDf['transcript_2']]))
        leftovers = unqXscriptSet - minXscripts
        oopsSet = set(pd.concat([oops1['transcript_1'], oops2['transcript_2']]))
        
        checkRemoved = leftovers == oopsSet
        
        if checkRemoved:
                print ("Valid Merge")
        else:
                print ("Missing transcripts in the final union")
        
        
        ergMatch = unionDf['ERG_id_' + args.GTF1] == unionDf['ERG_id_' + args.GTF2]
        unionDf['flag_ERG_match'] = ergMatch * 1
        
        unionDf['flag_ERG_noIR'] = ((unionDf['flag_ERG_match'] == 1) & (unionDf['flag_IR'] == 0)) * 1 
        unionDf['flag_ERG_wIR'] = ((unionDf['flag_ERG_match'] == 1) & (unionDf['flag_IR'] == 1)) * 1 

        unionDf = unionDf.drop(columns="merge_check")
        
        if args.colPrefix:
                colPrefix = args.colPrefix
                unionDf = unionDf.add_prefix(colPrefix + "_")
                
                unionDf = unionDf.rename(columns={colPrefix + "_transcript_1":"transcript_1",
                                colPrefix + "_transcript_2":"transcript_2"})
        
        
        outputFile = args.outDir + "/" + args.GTF1 + "_" + args.GTF2 + "_" + "map.csv"
        
        try:
                unionDf.to_csv(outputFile,index=False)
        except OSError:
                raise OSError("Output directory must already exist.")

        print ("Complete!")                
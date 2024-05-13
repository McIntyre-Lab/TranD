#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 15:05:58 2022

@author: nkeil
"""

import argparse
import pandas as pd
import numpy as np
import sys
import csv
import time
import trand.io

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Subset Reference GTF based on a list of transcript IDs")

    # Input arguments
    parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="GTF file where transcript_id and gene_id are the first two attributes in the attributes column (order does not matter)")
    parser.add_argument("-t", "--list-type", dest="inType", required=True, choices=['transcript_id','gene_id','chr'], default='transcript_id', help="Type of feature in inclusion/exclusion list provided (chr, transcript_id, gene_id, or chr, default is transcript_id)")
    parser.add_argument("-i", "--include-list", dest="inInclude", required=False, help="List of transcript or gene IDs to include in output from input GTF (no header) - cannot be used with exclusion list")
    parser.add_argument("-e", "--exclude-list", dest="inExclude", required=False, help="List of transcript or gene IDs to exclude in output from input GTF (no header) - cannot be used with inclusion list")
  
    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file name for subset GTF")

    args = parser.parse_args()
    return args


def main():
    # Get GTF file and count total lines

    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    # inInclude = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/design_files/df_noHeader_dmel6_sex_det_gene_lst_01ksb.csv"
    # listVar = 'gene_id'
    # outFile = ""
    # inExclude = None
    
    
    inGTF = args.inGTF
    
    # Get list ID variable
    listVar = args.inType

    inExclude = args.inExclude
    inInclude = args.inInclude

    outFile = args.outFile
    
    # Check that only one list argument was provided
    if inInclude is None and inExclude is None:
        print("ERROR: No list provided to subset GTF on")
        sys.exit()
    elif inInclude is not None and inExclude is None:
        include = True
        idDF = pd.read_csv(inInclude, names=[listVar])
    elif inInclude is None and inExclude is not None:
        include = False
        idDF = pd.read_csv(inExclude, names=[listVar])
    else:
        print("ERROR: Cannot use both inclusion and exclusion list - only provide one or the other")
        sys.exit()
    
  
    gtf = trand.io.read_exon_data_from_file(inGTF)    

    if include:
        # If list is for inclusion
        # Extract GTF elements that have a value contained in list of IDs
        filterDF = gtf[gtf[listVar].isin(idDF[listVar])]
    else:
        # If list is for exclusion
        # Extract GTF elements that do NOT have a value contained in list of IDs
        filterDF = gtf[~gtf[listVar].isin(idDF[listVar])]
    
    # Output filtered file
    filterDF[['chr','source','feature','start','end','score','strand','frame',
             'attribute']].to_csv(outFile,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)
                              
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()


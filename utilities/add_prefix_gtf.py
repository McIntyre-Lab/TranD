#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 13:14:17 2022

@author: nkeil
"""
import argparse
import pandas as pd
import numpy as np
import csv

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Add prefix to transcript id or gene id in gtf")

    # Input arguments
    parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="GTF file where transcript_id and gene_id are the first two attributes in the attributes column (order does not matter)")
    parser.add_argument("-t", "--id-type", dest="idType", required=True, choices=['transcript_id','gene_id'], default='transcript_id', help="Type of feature to which the prefix is added (transcript_id or gene_id, default is transcript_id)")
    parser.add_argument("-p", "--prefix", dest="inPrefix", required=True, help="Prefix to be added to transcript_id or gene_id")
    
    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file name for subset GTF")

    args = parser.parse_args()
    return args


##Get position of gene_id or transcript_id in list
def getIndex(gtf_id, attr_list):
    if gtf_id == "gene_id":
        for x in attr_list:
            if "gene_id" in x:
                return attr_list.index(x)
    elif gtf_id == "transcript_id":
        for x in attr_list:
            if 'transcript_id' in x:
                return attr_list.index(x)   

def main():
    # Get GTF file and count total lines
    gtf = pd.read_csv(args.inGTF,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], dtype=str, sep="\t",low_memory=False)
    
    
    # Get id type variable
    idVar = args.idType
    
    # Get prefix variable
    prefix = args.inPrefix
    
   
    
    # Get gene_id and transcript_id values from attributes column
    for i in gtf.index:
        raw_attrs = gtf.at[i, 'attribute']
        attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
        g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = np.nan, np.nan
        
        
                 
            
        for item in g_t_attrs:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
            elif 'transcript_id' in item:
                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
        if gene_id == np.nan:
            print("WARNING: gene_id not found in {}".format(gtf[i]))
        if transcript_id == np.nan and gtf.at[i, 'feature'] != "gene" :
            print("WARNING: transcript_id not found in {}".format(gtf[i]))
        # logger.debug("Gene: {}, Transcript: {}", gene_id, transcript_id)
        gtf.at[i, 'gene_id'] = str(gene_id)
        gtf.at[i, 'transcript_id'] = str(transcript_id)
        
        if idVar == "gene_id":
            prefix_gene_id = "gene_id " + '\"' + str(prefix) + "_" + str(gtf.at[i, 'gene_id']) + '\"'
        elif idVar == "transcript_id":
            prefix_transcript_id = "transcript_id " + '\"' + str(prefix) + "_" + gtf.at[i, 'transcript_id'] + '\"'
        
        prefix_attr_list=attr_list
        
        if idVar == "gene_id":
            prefix_attr_list[getIndex(idVar, attr_list)]=prefix_gene_id
            
        elif idVar == "transcript_id":
            if gtf.at[i, 'feature'] != "gene" :
                prefix_attr_list[getIndex(idVar, attr_list)]=prefix_transcript_id
        
        prefix_atrribute= ";".join(prefix_attr_list)
        
        gtf.at[i,'prefix_attribute']=prefix_atrribute
        
        
        
    # Output file with prefix
    gtf[['chr','source','feature','start','end','score','strand','frame',
             'prefix_attribute']].to_csv(args.outFile,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

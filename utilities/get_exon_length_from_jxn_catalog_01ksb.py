#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 17:24:49 2024

@author: k.bankole
"""

import argparse
import pandas as pd
import trand.io

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Gets the lengths of exons from TranD's junction catalog "
                                         "for every exon except the first and last (not possible using only junction information). "
                                         "Outputs a CSV that lists all exons for all transcripts and their lengths. "
                                         "If exon-length is 0, it is the first/last exon."
                                         "Add -t to change the nt threshold for mini-exons (number of mini-exons "
                                         "is output by print statement).")
        
        # INPUT
        parser.add_argument(
                "-jc",
                "--junction-catalog",
                dest="catalogFile",
                required=True,
                help="Location of \'junction_catalog.csv\' file from TranD output.")
        
        
        parser.add_argument(
                "-t",
                "--threshold",
                dest="miniThresh",
                default=8,
                required=False,
                help="Optional mini-exon nt size threshold. Default: 8")
        
        
        # OUTPUT
        parser.add_argument("-o",
                            "--output",
                            dest="outFile", 
                            required=True,
                            help="Path and filename for output CSV")
        
        args = parser.parse_args()
        return args

def main():
    
    # catalogFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dmel6_ujc_junction_catalog.csv"
    # miniThresh = 8
    
    catalogFile = args.catalogFile
    outFile = args.outFile
    miniThresh = args.miniThresh
    
    catDf = pd.read_csv(catalogFile, low_memory=False)    
    catDf = catDf.sort_values(by=['gene_id','transcript_id','coords'])
    
    xscriptDf = catDf.groupby(['gene_id','transcript_id']).agg(list).reset_index()
    
    coordDf = xscriptDf.explode('coords')
    
    coordDf[['chr', 'start', 'end', 'strand']] = coordDf['coords'].str.split(':', expand=True)
    
    coordDf['start'] = coordDf['start'].astype(int)
    coordDf['end'] = coordDf['start'].astype(int)
    
    coordDf = coordDf[['gene_id','transcript_id','start','end']]
    
    coordDf = coordDf.rename(columns={
        'gene_id':'geneID',
        'transcript_id':'transcriptID',
        })
    
    jxnDf = coordDf.groupby(['geneID','transcriptID']).apply(lambda x: [(start, end) for start, end in zip(x['start'], x['end'])]).reset_index(name='jxns')
    
    def process_jxns(jxns):
        # Extract start and end values from sorted list of tuples
        startValues, endValues = zip(*sorted(jxns))
        
        # Create list of exon tuples by skipping first start and last end
        exons = list(zip(endValues[:-1], startValues[1:]))
        
        # Define first and last exon tuples
        firstExon = ('firstDonor', startValues[0])
        lastExon = (endValues[-1], 'lastAcceptor')
        
        # Prepend first exon and append last exon to the exon list
        exons = [firstExon] + exons + [lastExon]
        return exons

    # Apply the process_jxns function to each row in the 'jxns' column
    jxnDf['exons'] = jxnDf['jxns'].apply(process_jxns)
    
    jxnDf = jxnDf[['geneID','transcriptID','exons']]
    
    jxnDf = jxnDf.explode('exons')
    outDf = jxnDf.copy()
    
    outDf[['start', 'end']] = pd.DataFrame(outDf['exons'].tolist(), index=jxnDf.index)
    outDf = outDf[['geneID','transcriptID','start','end']]
    
    outDf['exonLength'] = outDf.apply(lambda row: row['end'] - row['start'] if isinstance(row['end'], int) and isinstance(row['start'], int) else 0, axis=1)
    
    miniExons = outDf[(outDf['exonLength'] <= 8) & (outDf['exonLength'] > 0)]
    
    print("There are {} mini-exons in this output. (nt size <= {}.)".format(len(miniExons),miniThresh))
    outDf.to_csv(outFile,index=False)
    

if __name__ == '__main__':
        global args
        args = getOptions()
        main()
    
    
    
    
    
    
    
    
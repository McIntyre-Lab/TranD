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
        parser = argparse.ArgumentParser(description="Uses junction catalog to get lengths of exons."
                                         "If exon-length is 0, it is the first/last exon (not included in "
                                         "the junction catalogue.")
        
        # INPUT
        parser.add_argument(
                "-jc",
                "--junction-catalog",
                dest="catFile",
                required=True,
                help="Location of event analysis ER file from gene mode output"
                )
        
        
        # OUTPUT
        parser.add_argument("-o",
                            "--output",
                            dest="outFile", 
                            required=True,
                            help="Path and filename for output CSV")
        
        args = parser.parse_args()
        return args

def main():
    
    # catFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dmel6_ujc_junction_catalog.csv"
    
    catFile = args.catFile
    outFile = args.outFile
    
    
    jxnDf = pd.read_csv(catFile, low_memory=False)

    # test = jxnDf[jxnDf['gene_id'] == 'FBgn0000504'].copy()
    
    jxnDf = jxnDf.sort_values(by=['gene_id','transcript_id','coords'])
    
    grp = jxnDf.groupby(['gene_id','transcript_id']).agg(list).reset_index()
    
    newExonDf = grp.explode('coords')
    
    newExonDf[['chr', 'start', 'end', 'strand']] = newExonDf['coords'].str.split(':', expand=True)
    
    newExonDf['start'] = newExonDf['start'].astype(int)
    newExonDf['end'] = newExonDf['start'].astype(int)
    
    newExonDf = newExonDf[['gene_id','transcript_id','start','end']]
    
    newExonDf = newExonDf.rename(columns={
        'gene_id':'geneID',
        'transcript_id':'transcriptID',
        })
    
    result = newExonDf.groupby(['geneID','transcriptID']).apply(lambda x: [(start, end) for start, end in zip(x['start'], x['end'])]).reset_index(name='jxns')
    
    def process_jxns(jxns):
        # Extract start and end values from sorted list of tuples
        startValues, endValues = zip(*sorted(jxns))
        
        # Create list of exon tuples by skipping first start and last end
        exons = list(zip(endValues[:-1], startValues[1:]))
        
        # Define first and last exon tuples
        firstExon = ('firstDonor', startValues[0])
        lastExon = (endValues[-1], 'lastAcceptor')
        
        # firstExon = ('0', startValues[0])
        # lastExon = (endValues[-1], '0')
        
        # Prepend first exon and append last exon to the exon list
        exons = [firstExon] + exons + [lastExon]
        
        return exons

    # Apply the process_jxns function to each row in the 'jxns' column
    result['exons'] = result['jxns'].apply(process_jxns)
    
    result = result[['geneID','transcriptID','exons']]
    
    result = result.explode('exons')
    newNewDf = result.copy()
    
    newNewDf[['start', 'end']] = pd.DataFrame(newNewDf['exons'].tolist(), index=result.index)
    newNewDf = newNewDf[['geneID','transcriptID','start','end']]
    
    newNewDf['exonLength'] = newNewDf.apply(lambda row: row['end'] - row['start'] if isinstance(row['end'], int) and isinstance(row['start'], int) else 0, axis=1)
    
    # miniExons = newNewDf[(newNewDf['exonLength'] <= 8) & (newNewDf['exonLength'] > 0)]
    
    newNewDf.to_csv(outFile,index=False)
    

if __name__ == '__main__':
        global args
        args = getOptions()
        main()
    
    
    
    
    
    
    
    
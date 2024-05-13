#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io
import re
import time

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Uses \"event_analysis_er.csv\" from TranD 1 GTF gene mode"
                                         "to, per gene, flag what exon regions a transcript has. "
                                         "Outputs a CSV in stacked form (one row per "
                                         "exon region per transcript.")
        
        # INPUT
        parser.add_argument(
                "-ea",
                "--event-analysis",
                dest="eaFile",
                required=True,
                help="Location of event analysis ER file from gene mode output"
                )
        
        # OUTPUT
        parser.add_argument("-o",
                            "--output-csv",
                            dest="outFile", 
                            required=True,
                            help="Output file path")
        
        args = parser.parse_args()
        return args


def main():
    
    # PROJ="/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/test_gene_based_ERGs"
    PROJ="/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs"
    
    # eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_ea_LOC110191367_noIR.csv"
    # eaFile = "{}/sub_keepIR_event_analysis_er.csv".format(PROJ)
    # eaFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dser1_ujc_event_analysis_er.csv"
    # outFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/test_ER_flags/test_ER_flags.csv"
    
    eaFile = args.eaFile
    outFile = args.outFile

    eaDf = pd.read_csv(eaFile, low_memory=False)       
    
    infoDf = eaDf[['er_id','er_transcript_ids','gene_id']].copy()
    
    infoDf.columns = ['ER','jxnHash','geneID']
    infoDf['jxnHash'] = infoDf['jxnHash'].str.split('|')
    infoDf = infoDf.explode('jxnHash')
    
    geneDct = dict(infoDf.groupby('geneID').apply(lambda x: set(x['ER'])))
    xscriptDct = dict(infoDf.groupby('jxnHash').apply(lambda x: set(x['ER'])))
    
    loopLst = [tuple(x) for x in infoDf[['geneID','jxnHash']].drop_duplicates().to_records(index=False)]
        
    xscriptLst = []
    geneLst = []
    erLst = []
    flagLst = []
    for gene, transcript in loopLst:
        
        # gene = row['geneID']
        # transcript = row['jxnHash']
        
        geneERSet = geneDct.get(gene)
        
        geneERSet = sorted(geneERSet, key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("exon_")[1]))
        
        xscriptERSet = xscriptDct.get(transcript)
        
        for exonRegion in geneERSet:
            
            if exonRegion in xscriptERSet:
                flag = 1 
            else:
                flag = 0
            
            xscriptLst.append(transcript)
            geneLst.append(gene)
            erLst.append(exonRegion)
            flagLst.append(flag)
    
    outDf = pd.DataFrame({
        'jxnHash':xscriptLst,
        'geneID':geneLst,
        'exonRegion':erLst,
        'flag_ERPresent':flagLst
    })
    
    outDf.to_csv(outFile,index=False)
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

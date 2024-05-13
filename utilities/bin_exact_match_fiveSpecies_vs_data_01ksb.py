#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 12:08:24 2024

@author: k.bankole
"""

import argparse
import pandas as pd

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Take in data uniq jxnHash CSV and corresponding "
                                         "fiveSpecies flag file and identify shared jxnHashes between "
                                         "the data and the reference. Output a list of jxnHashes that "
                                         "are exact matches.")
        
        # INPUT
        parser.add_argument(
                "-f",
                "--flag-file",
                dest="inFlag",
                required=True,
                help="Location of fiveSpecies flag file"
        )
        parser.add_argument(
                "-d",
                "--data-csv",
                dest="inData",
                required=True,
                help="Location of data unique jxnHash CSV file"
        )
        
        
        
        # OUTPUT
        parser.add_argument("-o",
                            "--output",
                            dest="outFile", 
                            required=True,
                            help="Output list path")
        
        args = parser.parse_args()
        return args


def main():
    
    # inFlag = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/flag_fiveSpecies_2_dmel6_ujc.csv"
    # inData = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/mel_2_dmel6_uniq_jxnHash.csv"
    # outFile = ""
    
    inFlag = args.inFlag
    inData = args.inData
    outFile = args.outFile
    
    flagDf = pd.read_csv(inFlag,low_memory=False)    
    dataDf = pd.read_csv(inData,low_memory=False)    


    dataDf = dataDf[['jxnHash']]
    
    flagDf = flagDf[[col for col in flagDf.columns if 'jxnHash' in col]]
    flagDf.columns = ['jxnHash']
    
    
    refJxnHashLst = flagDf['jxnHash'].unique()
    
    exactMatchDf = dataDf[dataDf['jxnHash'].isin(refJxnHashLst)]
    
    exactMatchDf.to_csv(outFile,index=False,header=False)
    
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

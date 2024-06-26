#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")

    # parser.add_argument(
    #     "-p",
    #     "--pattern-file",
    #     dest="erpFile",
    #     required=True,
    #     help="Location of ER pattern file"
    # )

    # parser.add_argument(
    #     "-c",
    #     "--count-file",
    #     dest="countFile",
    #     required=False,
    #     help="Location of counts per jxnHash"
    # )

    # # # Output data
    # parser.add_argument(
    #     "-o",
    #     "--outdir",
    #     dest="outdir",
    #     required=True,
    #     help="Output directory"
    # )

    # parser.add_argument(
    #     "-n",
    #     "--data-filename",
    #     dest="fileName",
    #     required=True,
    #     help="Name of data GTF for output files. Required."
    # )

    # parser.add_argument(
    #     "-x",
    #     "--prefix",
    #     dest="prefix",
    #     required=False,
    #     help="Prefix for output files."
    # )

    args = parser.parse_args()
    return args

def main():
    
    inAnnoERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er_vs_fiveSpecies_2_dmel6_ujc_ERP.csv"
    inDataERPFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/compare_fiveSpecies_er_vs_data_gtf/fiveSpecies_2_dmel6_ujc_er_vs_mel_2_dmel6_uniq_jxnHash_ERP.csv"
    in5SpFlagFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/flag_fiveSpecies_2_dmel6_ujc.csv"
    
    inAnnoERPDf = pd.read_csv(inAnnoERPFile, low_memory=False) 
    inDataERPDf = pd.read_csv(inDataERPFile, low_memory=False) 
    in5SpFlagDf = pd.read_csv(in5SpFlagFile, low_memory=False) 
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    # main()

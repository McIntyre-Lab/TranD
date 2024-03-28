#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 12:31:29 2024

@author: k.bankole
"""

import argparse
import pandas as pd
import trand.io


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="")
    
    # Input data
    parser.add_argument("-ud",
                        "--ujc-dsc",
                        dest="ujcDsc", 
                        required=True,
                        help="Path to catted UJC dsc file (created using cat_and_count_ujc_dsc_01ksb.py)")
        
    # Output data
    parser.add_argument("-o",
                        "--output-gtf",
                        dest="outGTF", 
                        required=True,
                        help="Path and filename for output GTF")
    
    
    args = parser.parse_args()
    return args

def main():
    
    # # dscFile = "/nfshome/k.bankole/Desktop/test_dsc/out_ujc_dsc.csv"
    # dscFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/rmg_lmm_dros_data/test.csv"
    # gtfOutPath = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_ujc_out_to_gtf/test.gtf"

    dscFile = args.ujcDsc
    gtfOutPath = args.outGTF
    
    dscDf = pd.read_csv(dscFile, low_memory=False)
    
    infoDf = dscDf[['jxnHash','chr','strand','donorStart','acceptorEnd','jxnString']].copy(deep=True)
    infoDf['chr'] = "chromosome" + infoDf['chr']
    infoDf['jxnString'] = "chromosome" + infoDf['jxnString']
    
    print ("Number of input jxnHash: {}".format(len(infoDf['jxnHash'])))
    print ("Number of input unique jxnHash: {}".format(infoDf['jxnHash'].nunique()))

    seqnameLst = []
    startLst = []
    endLst = []
    hashLst = []
    strandLst = []
    geneIDLst = []    
    
    for row in infoDf.to_dict('records'):
            seqname = row['chr']
            strand = row['strand']
            jxnHash = row['jxnHash']
            geneID = row['jxnHash']
            
            firstStart = row['donorStart']
            lastEnd = row['acceptorEnd']
            
            jxnString = row['jxnString']
            flagMono = 'monoexon' in row['jxnString']
            
            if flagMono:
                seqnameLst.append(seqname)
                startLst.append(firstStart)
                endLst.append(lastEnd)
                hashLst.append(jxnHash)
                strandLst.append(strand)
                geneIDLst.append(geneID)
                
            else:
                
                jxnLst = jxnString.split(seqname + '_')[1].split(strand + '_')[1].split('_')
                
                seqnameLst.append(seqname)
                startLst.append(firstStart)
                endLst.append(jxnLst[0])
                hashLst.append(jxnHash)
                strandLst.append(strand)
                geneIDLst.append(geneID)
                
                for i in range(1, len(jxnLst) -1):                
                    seqnameLst.append(seqname)
                    startLst.append(jxnLst[i])
                    endLst.append(jxnLst[i+1])
                    hashLst.append(jxnHash)
                    strandLst.append(strand)
                    geneIDLst.append(geneID)
                
                seqnameLst.append(seqname)
                startLst.append(jxnLst[-1])
                endLst.append(lastEnd)
                hashLst.append(jxnHash)
                strandLst.append(strand)
                geneIDLst.append(geneID)

    
    outExonDf = pd.DataFrame(
            {
                    'seqname':seqnameLst,
                    'start':startLst,
                    'end':endLst,
                    'strand':strandLst,
                    'transcript_id':hashLst,
                    'gene_id':geneIDLst
            })
    
    numColumns = ['start','end']
    outExonDf[numColumns] = outExonDf[numColumns].astype(int)
    outExonDf = outExonDf.sort_values(by=['seqname','transcript_id','start'])
    result = outExonDf[outExonDf['end'] < outExonDf['start']]
    
    print ("Number of output unique jxnHash: {}".format(outExonDf['transcript_id'].nunique()))
    
    trand.io.write_gtf(data=outExonDf, out_fhs={"gtf":gtfOutPath}, fh_name="gtf")
    
    
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
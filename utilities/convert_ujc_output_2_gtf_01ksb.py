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
    # parser.add_argument("-g",
    #                     "--input-gtf",
    #                     dest="inGTF", 
    #                     required=True,
    #                     help="GTF to find genes with one transcript")
    
    # Output data
    parser.add_argument("-", "--", dest="", required=True, help="")
    
    args = parser.parse_args()
    return args

def main():
    
    
    indexFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_ujc_out_to_gtf/head50_dmel650_2_dmel6_ujc_xscript_index.csv"
    idFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_ujc_out_to_gtf/head50_dmel650_2_dmel6_ujc_id.csv"
    gtfOutPath = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_conv_ujc_out_to_gtf/test.gtf"

    indexDf = pd.read_csv(indexFile, low_memory=False)
    idDf = pd.read_csv(idFile, low_memory=False)
    
    groupDf = indexDf.groupby('jxnHash').agg({'jxnString':set})
    
    # Check collisions
    a = groupDf[groupDf['jxnString'].apply(lambda x: len(x) > 1)]
    if len(a) > 1:
        print ("Collision!!!")
    
    
    groupDf['jxnString'] = groupDf['jxnString'].apply(lambda x: str(list(x)[0]))
    # Need jxnHash, jxnString, and start end
    # Helpful: chr, strand, source
    
    mergeDf = pd.merge(idDf, groupDf, on='jxnHash', how='outer',indicator='merge_check')
    
    infoDf = mergeDf[['jxnHash','chr','strand','donorStart','acceptorEnd','jxnString']].copy(deep=True)
    
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
                
                jxnLst = jxnString.split("_")[2:]
                
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
    
    outExonDf = outExonDf.sort_values(by=['start','seqname'])
    
    print ("Number of output unique jxnHash: {}".format(outExonDf['transcript_id'].nunique()))
    
    
    # inGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    
    
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
#!/usr/bin/env python

import argparse
import pandas as pd

import trand.io
import os

# def getOptions():
#     # Parse command line arguments
#     parser = argparse.ArgumentParser(description="")

#     # Input data
#     parser.add_argument("-", "--", dest="", required=True, help="")

#     # Output data
#     parser.add_argument("-", "--", dest="", required=True, help="")

#     args = parser.parse_args()
#     return args

def main():    

    pairFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_dmel6_dsim2_dser1_sexDet/dser1_pairwise_transcript_distance.csv"
    inDf = pd.read_csv(pairFile, low_memory=False)
    
    pairDf1 = inDf[inDf['transcript_1'].str.contains('307278e675409195678d40133e74781215e4950ca58fe890b0cb2b2dffb8d652')].copy()
    pairDf2 = inDf[inDf['transcript_2'].str.contains('307278e675409195678d40133e74781215e4950ca58fe890b0cb2b2dffb8d652')].copy()
    
    pairDf1 = pairDf1[['transcript_1','transcript_2','ER_only_T1','ER_ovlp']].fillna('')
    pairDf1.columns = ['transcriptID','otherTranscript','ER_only','ER_ovlp']
    
    pairDf2 = pairDf2[['transcript_2','transcript_1','ER_only_T2','ER_ovlp']].fillna('')
    pairDf2.columns = ['transcriptID','otherTranscript','ER_only','ER_ovlp']
     
    pairDf = pd.concat([pairDf1,pairDf2]).reset_index(drop=True)
    
    pairDf['ERs'] = pairDf['ER_only'] + '|' + pairDf['ER_ovlp']
    
    # eaFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dser1_ujc_event_analysis_er.csv"
    # inDf = pd.read_csv(eaFile, low_memory=False)

    # eaDf = inDf[inDf['gene_id'].str.contains('LOC110191367')].copy()
    
    seqnameLst = []
    startLst = []
    endLst = []
    strandLst = []
    ergIDLst = []
    geneIDLst = []
    
    rowNum = 0
    for row in pairDf.to_dict('records'):
        rowNum += 1
        
        erList = row['ERs'].split('|')
        
        
        for er in erList:
            if er != '':
                    seqname = er.split(':')[0]
                    start = int(er.split(':')[1])
                    end = int(er.split(':')[2])
                    strand = er.split(':')[3]
                    
                    seqnameLst.append(seqname)
                    strandLst.append(strand)
                    ergIDLst.append(f'row_{rowNum}')
                    geneIDLst.append(row['otherTranscript'])
                    
                    startLst.append(start)
                    endLst.append(end)
                
                
                
    exonDf = pd.DataFrame(
            {
                    'seqname':seqnameLst,
                    'start':startLst,
                    'end':endLst,
                    'strand':strandLst,
                    'transcript_id':ergIDLst,
                    'gene_id':geneIDLst
            })

    exonDf = exonDf.sort_values(by=['seqname','transcript_id','start']).reset_index(drop=True)

    gtf_output_file = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/out.gtf"
    if os.path.isfile(gtf_output_file):
            os.remove(gtf_output_file)
    trand.io.write_gtf(data=exonDf, out_fhs={"gtf":gtf_output_file}, fh_name="gtf")



if __name__ == '__main__':
    # Parse command line arguments
    global args
    # args = getOptions()
    main()

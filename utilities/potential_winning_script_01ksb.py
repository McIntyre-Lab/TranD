#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:29:04 2024

@author: k.bankole
"""

import argparse
import pandas as pd
import trand.io
import time
import numpy as np

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Placeholder")
        
        # INPUT
        parser.add_argument(
                "-er",
                "--er-gtf",
                dest="erFile",
                required=True,
                help="Location of ER GTF"
                )
        
        
        parser.add_argument(
                "-d",
                "--data-gtf",
                dest="dataFile",
                required=True,
                help="Location of data GTF"
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
    
    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_data.gtf"    
    
    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/dmel_sexdet_fiveSpecies_vs_data/mel_2_dmel6_uniq_jxnHash_sexDet.gtf"
    
    outFile = "/nfshome/k.bankole/Desktop/test_folder"
    
    erFile = args.erFile
    dataFile = args.dataFile
    outFile = args.outFile
    
    alphatic = time.perf_counter()
    
    geneDf = trand.io.read_exon_data_from_file(erFile)
    dataDf = trand.io.read_exon_data_from_file(dataFile)
    
    geneDf = geneDf[['gene_id','seqname','start','end','strand']]
    geneDf = geneDf.sort_values(['seqname','gene_id','start'])
    
    geneDf['ER'] = geneDf['gene_id'] + ':ER' + (geneDf.groupby('gene_id').cumcount() + 1).astype(str)
    geneDct = dict(geneDf.groupby('gene_id').apply(lambda x: sorted(set(x['ER']), key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("exon_")[1]))))
    
    # TODO: CHECK THAT ALL SETS ARE OF SIZE ONE
    erDct = geneDf.groupby('ER').agg('first').to_dict(orient='index')
    
    # row = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+1])
    # row2 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+2])
    # row3 = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+3])
    # yourBoat = pd.DataFrame({'gene_id':'test','seqname':'test','start':0,'end':0,'strand':13341}, index=[len(dataDf)+4])
    
    # dataDf = pd.concat([dataDf,row,row2,row3,yourBoat])
    
    dataOnlyGnLst = []
    
    dataDf['numExon'] = dataDf.groupby('transcript_id')['transcript_id'].transform('count')

    records = dataDf.to_dict('records')
    
    for row in records:
        gene = row['gene_id']
        # jxnHash = row['transcript_id']
        
        matchingERIDLst = []
        
        if gene in geneDct.keys():
            for erID in geneDct.get(gene):
                # print(erID)
                erInfo = erDct.get(erID)
                # print(erInfo)
                
                # print("looping...")
                
                if max(row['start'], erInfo['start']) < min(row['end'], erInfo['end']):
                    # print ("Yippee!")
                    # print(row)
                    # print(erID)
                    # print(erInfo)
                    matchingERIDLst.append(erID)
                    
                    
            if matchingERIDLst:
                row['ER'] = matchingERIDLst
            else:
                row['dataOnlyER'] = "{}:{}_{}".format(gene,row['start'], row['end'])
        else:
            # TODO: What if the gene only appears in the reference? do we care?
            dataOnlyGnLst.append(row)
    
    dataWithERDf = pd.DataFrame(records)
    
    # TODO: May not be true but makes sense right? if there is an exon that overlaps two exon regions
    # it may not be true biological IR but its definitely overlapping an intron...
    dataWithERDf['flagIR'] = dataWithERDf['ER'].apply(lambda x: x if not type(x) is list else 1 if len(x) > 1 else 0)
    dataWithERDf['IRERs'] = dataWithERDf.apply(lambda x: tuple(x['ER']) if x['flagIR'] == 1 else np.nan, axis=1) 

    intmdDf = dataWithERDf[['transcript_id','ER','dataOnlyER','flagIR','IRERs','numExon']]
    intmdDf = intmdDf.explode('ER')
    
    xscriptERDf = intmdDf.groupby('transcript_id').agg({
        'ER':lambda x: set(x.dropna()),
        'numExon':max,
        'flagIR':max,
        'dataOnlyER':lambda x: set(x.dropna()),
        'IRERs':lambda x: set(tuple(sum(x.dropna(),())))
    }).reset_index()
    
    # Accounts for situations where transcripts have multiple IR events
    xscriptERDf['numIREvent'] = xscriptERDf.groupby('transcript_id')['flagIR'].transform('sum')
    
    xscriptERDct = dict(zip(xscriptERDf['transcript_id'], xscriptERDf['ER']))
    
    loopLst = [tuple(x) for x in dataWithERDf[['gene_id','transcript_id']].drop_duplicates().to_records(index=False)]
    
    xscriptLst = []
    geneLst = []
    erLst = []
    flagLst = []
    
    binaryDct = dict()
    for gene, transcript in loopLst:
        
        # gene = row['geneID']
        # transcript = row['jxnHash']
        
        geneERLst = geneDct.get(gene)
        xscriptERSet = xscriptERDct.get(transcript)
    
        binary = [1 if ER in xscriptERSet else 0 for ER in geneERLst]
        binary = ''.join(map(str, binary))
        
        binaryDct[transcript] = binary
        
        for exonRegion in geneERLst:
            
            if exonRegion in xscriptERSet:
                flag = 1 
            else:
                flag = 0
            
            xscriptLst.append(transcript)
            geneLst.append(gene)
            erLst.append(exonRegion)
            flagLst.append(flag)

    # Add Data Only ERs... somehow...
    outDf = pd.DataFrame({
        'jxnHash':xscriptLst,
        'geneID':geneLst,
        'exonRegion':erLst,
        'flag_ERPresent':flagLst
    })
    
    binaryDf = pd.DataFrame.from_dict(binaryDct, orient='index').reset_index()
    binaryDf.columns=['transcript_id','ERP']
    
    otherOutDf = pd.merge(xscriptERDf,binaryDf,on=['transcript_id'],how='outer',indicator='merge_check')
    
    otherOutDf['ERP'] = otherOutDf.apply(lambda x: x['ERP'] + str('1' * len(x['dataOnlyER'])) if x['dataOnlyER'] else x['ERP'],axis=1)
    
    otherOutDf['numER'] = otherOutDf['ER'].apply(len)
    otherOutDf['numDataOnlyER'] = otherOutDf['dataOnlyER'].apply(len)
    
    if not (otherOutDf['merge_check'] == 'both').all():
        print ("Something went wrong")
        # quit()
    
    otherOutDf = otherOutDf[['transcript_id','ERP','numExon','numER','numDataOnlyER','flagIR','numIREvent','IRERs']]
    
    otherOutDf['flagReverseIR'] = otherOutDf.apply(lambda x: 1 if x['numExon'] > x['numER'] + x['numDataOnlyER'] else 0, axis=1)
    
    otherOutDf['IRERs'] = otherOutDf['IRERs'].apply(lambda x: '|'.join(x) if x else np.nan)
    
    # Do not uncomment.
    # wideDf = pd.pivot_table(outDf, values='flag_ERPresent', index=['jxnHash','geneID'], columns='exonRegion', fill_value=0)
    
    outDf.to_csv(outFile,index=False)
    # Output data only gene list (remove from existing dataframe)
    
    omegatoc = time.perf_counter()
    
    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")        
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
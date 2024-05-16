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
    
    # erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
    # dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/FBgn0000662_data.gtf"    
    
    
    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/trand_2gtf_fiveSpecies_er_vs_data/fiveSpecies_2_dmel6_ujc_er.gtf"
    dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/dmel_sexdet_fiveSpecies_vs_data/mel_2_dmel6_uniq_jxnHash_sexDet.gtf"
    
    
    erFile = args.erFile
    dataFile = args.dataFile
    
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
    
    dataOnlyLst = []
    
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
                row['ER(s)'] = matchingERIDLst
            else:
                row['ER(s)'] = "{}:{}_{}".format(gene,row['start'], row['end'])
        else:
            dataOnlyLst.append(row)
                    
    
    dataWithERDf = pd.DataFrame(records).dropna(subset=['ER(s)'])
    
    
    # Output data only gene list (remove from existing dataframe)
    
    # What if the gene only appears in the reference? does that matter? 
    
    omegatoc = time.perf_counter()
    
    print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")    
    
    
    # May not be true but makes sense right? if there is an exon that overlaps two exon regions
    # it may not be true biological IR but its definitely overlapping an intron...
    
    # dataWithERDf['flag_IR'] = dataWithERDf.apply(1 if len(dataWithERDf['ER(s)']) > 1 else 0)
    
    # dataDf.groupby('gene_id')
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
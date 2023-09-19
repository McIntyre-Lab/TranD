#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:15:16 2023

@author: k.bankole
"""

import pandas as pd
import time
import numpy as np
from collections import Counter


if __name__ == '__main__':
        print ("Loading...")

        omegatic = time.perf_counter()
        
        longfile = "/nfshome/k.bankole/Desktop/1GTFtest/GCF_000090745.1_AnoCar2.0_genomic.gtf"
        shortfile = "/nfshome/k.bankole/Desktop/1GTFtest/testgtf.gtf"
        
        long = True
        
        if long:
                data = pd.read_csv(longfile, sep='\t', comment='#', header=None, low_memory=False)
        else:
                data = pd.read_csv(shortfile, sep='\t', comment='#', header=None, low_memory=False)

        
        all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                           'attributes', 'comments']
        
        drop_columns = ['source', 'feature', 'score', 'frame', 'comments']

        file_cols = data.columns

        if len(file_cols) < len(all_gtf_columns):
            gtf_cols = all_gtf_columns[:len(file_cols)]
        data.columns = gtf_cols
        drop_cols = [x for x in drop_columns if x in gtf_cols]
        
        data = data[data['feature'] == 'exon']
        data = data.drop(labels=drop_cols, axis=1)

        data['seqname'] = data['seqname'].astype(str)
        data['start'] = data['start'].astype(int)
        data['end'] = data['end'].astype(int)
        
        data.reset_index(drop=True, inplace=True)
        
        seqnameLst = []
        startLst = []
        endLst = []
        strandLst = []
        geneIDLst = []
        xscriptIDLst = []
        
        for row in data.to_dict('records'):
                rawAttr = row['attributes']
                attrLst = [x.strip() for x in rawAttr.strip().split(';')]
                gnTrAttr = [x for x in attrLst if 'transcript_id' in x or 'gene_id' in x]
                gene_id, transcript_id = None, None
                
                for item in gnTrAttr:
                        if 'gene_id' in item:
                                gene_id = item.split('gene_id')[1].strip().strip('\"')
                        elif 'transcript_id' in item:
                                transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
                                
                if not gene_id:
                        print("gene_id not found in '{}'", row)
                        gene_id = None
                        
                if not transcript_id:
                        print("transcript_id not found in '{}'", row)
                        transcript_id = None

                
                seqnameLst.append(row['seqname'])
                startLst.append(row['start'])
                endLst.append(row['end'])
                strandLst.append(row['strand'])
                
                geneIDLst.append(gene_id)
                xscriptIDLst.append(transcript_id)
                
        
        newData = pd.DataFrame(
                {
                        'seqname':seqnameLst,
                        'start':startLst,
                        'end':endLst,
                        'strand':strandLst,
                        'gene_id':geneIDLst,
                        'transcript_id':xscriptIDLst
                })
        
        a = newData.groupby('gene_id')['transcript_id']
        
        print("Exon data rows: {}".format(newData.shape[0]))
        
        missing_value_num = newData.isnull().sum().sum()
        if missing_value_num > 0:
                print("Total number of missing values: {}", missing_value_num)
        else:
                print("No missing values in data")
        
        gene_id_missing_value_num = newData['gene_id'].isnull().sum()
        
        transcript_id_missing_value_num = newData['transcript_id'].isnull().sum()
        
        if gene_id_missing_value_num > 0:
                print("Missing gene_id value number: {}", gene_id_missing_value_num)
        if transcript_id_missing_value_num > 0:
                print("Missing transcript_id value number: {}", transcript_id_missing_value_num)
        
                
        newData['start'] = pd.to_numeric(newData['start'], downcast="unsigned")
        newData['end'] = pd.to_numeric(newData['end'], downcast="unsigned")
        
        toc = time.perf_counter()       
        
        print(f"Complete,  took {toc-omegatic:0.4f} seconds.")
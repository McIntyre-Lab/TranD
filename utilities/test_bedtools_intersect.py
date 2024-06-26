#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 13:24:57 2024

@author: k.bankole
"""

import trand.io
import pandas as pd
import numpy as np
import time

from pybedtools import BedTool

import trand.bedtools as BD

erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_fiveSpecies_2_dmel6_ujc_er.gtf"
dataFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/FBgn0000662_data.gtf"


alphatic = time.perf_counter()

inGeneDf = trand.io.read_exon_data_from_file(erFile)
inDataDf = trand.io.read_exon_data_from_file(dataFile)

uniqDataGeneSet = set(inDataDf['gene_id'])
uniqRefGeneSet = set(inGeneDf['gene_id'])

refOnlyGnLst = list(uniqRefGeneSet - uniqDataGeneSet)
dataOnlyGnLst = list(uniqDataGeneSet - uniqRefGeneSet)

genesInBoth = list(uniqRefGeneSet.intersection(uniqDataGeneSet))

geneDf = inGeneDf[inGeneDf['gene_id'].isin(genesInBoth)].copy()
dataDf = inDataDf[inDataDf['gene_id'].isin(genesInBoth)].copy()

geneDf = geneDf[['gene_id', 'seqname', 'start', 'end', 'strand']].copy()
geneDf = geneDf.sort_values(['seqname', 'gene_id', 'start'])

geneDf['ER'] = geneDf['gene_id'] + ':ER' + \
    (geneDf.groupby('gene_id').cumcount() + 1).astype(str)

singleStrandGene = geneDf.groupby('gene_id').agg(
    set)['strand'].apply(lambda x: len(x) == 1)

if not singleStrandGene.all():
    print("There are genes belonging to more than one strand. Quitting.")
    quit()

geneDct = dict(geneDf.groupby('gene_id').apply(lambda x:
                                               sorted(set(x['ER']), key=lambda x: int(x.split("ER")[
                                                      1]) if 'ER' in x else int(x.split("exon_")[1]))
                                               if (x['strand'] == "+").all()
                                               else
                                               sorted(set(x['ER']), key=lambda x: int(x.split("ER")[
                                                   1]) if 'ER' in x else int(x.split("exon_")[1]), reverse=True)))

erDct = geneDf.groupby('ER').agg('first').to_dict(orient='index')





test = BD.prep_bed_for_ea(dataDf)












myWay()

    
def myWay():
    
    print("Looping")
    tic = time.perf_counter()
    
    dataDf['numExon'] = dataDf.groupby('transcript_id')[
        'transcript_id'].transform('count')
    
    dataDf['dataOnlyExon'] = np.nan

    
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
                    # print(row)
                    # print(erID)
                    # print(erInfo)
    
                    matchingERIDLst.append(erID)
    
            if matchingERIDLst:
                row['ER'] = matchingERIDLst
            else:
                row['dataOnlyExon'] = "{}:{}_{}".format(
                    gene, row['start'], row['end'])
    
    dataWithERDf = pd.DataFrame(records)
    
    toc = time.perf_counter()
    
    print(f"Loop Complete! Took {(toc-tic):0.4f} seconds.")




# omegatoc = time.perf_counter()

# print(f"Complete! Took {(omegatoc-alphatic):0.4f} seconds.")

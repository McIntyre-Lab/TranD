#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.

"""
Read and write files
"""

import csv
from numpy import nan
import pandas as pd
from loguru import logger
from pathlib import Path
import os
import time


def prepare_outdir(out_path, force):
    outdir = Path(out_path)
    logger.debug("Output directory: {}", str(outdir))
    if outdir.exists():
        if not force:
            exit("Not overwriting existing output directory without -f|--force. Exiting.")
    else:
        outdir.mkdir(parents=True, exist_ok=True)


def open_output_files(outdir, outfiles):
    """Open all necessary output files and return filehandles"""
    out_fhs = {}
    try:
        for outfile in outfiles:
            fh = Path(outdir) / outfiles[outfile]
            
            if outfile != 'og_fh':
                    open(fh, 'w')
            else:
                    if os.path.exists(fh):
                            os.remove(fh)
                            
                            
                    
            out_fhs[outfile] = fh
    except SystemError as e:
        logger.error("Cannot open one or more output file for writing: {}", e)
    return out_fhs


def write_output(data, out_fhs, fh_name):
    """Write results of event analysis to output files."""
    data.to_csv(out_fhs[fh_name], mode='a', header=False, index=False)


def get_gtf_attribute(transcript_id, gene_id):
    return f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'

def write_txt(txtLst, out_fhs, fh_name):
        """Write two column text output """
        with out_fhs[fh_name].open("a") as f:
                for line in txtLst:
                        f.write("{:<20} {:<20}\n".format(str(line[0]), str(line[1])))
                        # f.write((txtLst[0] + "{:<5} *2 " txtlst[1]).format(*row))
        
def write_gtf(data, out_fhs, fh_name):
    """Write output gtf files."""
    if len(data) == 0:
        return
    else:
        data.loc[:, 'source'] = "TranD"
        data.loc[:, 'feature'] = "exon"
        data.loc[:, 'score'] = "."
        data.loc[:, 'frame'] = "."
        data.loc[:, 'attribute'] = data.apply(lambda x: get_gtf_attribute(x['transcript_id'],
                                              x['gene_id']), axis=1)
        output_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                               'frame', 'attribute']
        data = data.reindex(columns=output_column_names)
        data.to_csv(out_fhs[fh_name], sep="\t", mode='a', index=False, header=False,
                    doublequote=False, quoting=csv.QUOTE_NONE)


# def read_exon_data_from_file(infile):
#     """
#     Create a pandas dataframe with exon records from a gtf file
#     Raw gtf data:
#     seqname source  feature  start end  score strand frame attributes
#     2L FlyBase 5UTR  7529  7679 . +  . gene_symbol "CG11023"; transcript_id "FBtr...
#     2L FlyBase exon  7529  8116 . +  . gene_symbol "CG11023"; transcript_id "FBtr...

#     Exon Fragment Data for Event Analysis:
#     seqname start end   strand gene_id      transcript_id
#     2L      7529  8116  +      FBgn0031208  FBtr0300689
#     2L      8193  9484  +      FBgn0031208  FBtr0300689
#     2L      9839  11344 -      FBgn0002121  FBtr0078169
#     """
#     all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
#                        'attributes', 'comments']
#     drop_columns = ['source', 'feature', 'score', 'frame', 'comments']
#     data = pd.read_csv(infile, sep='\t', comment='#', header=None, low_memory=False)
#     file_cols = data.columns
#     if len(file_cols) < len(all_gtf_columns):
#         gtf_cols = all_gtf_columns[:len(file_cols)]
#     data.columns = gtf_cols
#     drop_cols = [x for x in drop_columns if x in gtf_cols]
#     data['seqname'] = data['seqname'].astype(str)
#     data['start'] = data['start'].astype(int)
#     data['end'] = data['end'].astype(int)
#     data = data[data['feature'] == 'exon']
#     data = data.drop(labels=drop_cols, axis=1)
#     data.reset_index(drop=True, inplace=True)
#     for i in data.index:
#         raw_attrs = data.at[i, "attributes"]
#         attr_list = [x.strip() for x in raw_attrs.strip().split(';')]
#         g_t_attrs = [x for x in attr_list if 'transcript_id' in x or 'gene_id' in x]
#         gene_id, transcript_id = nan, nan
#         for item in g_t_attrs:
#             if 'gene_id' in item:
#                 gene_id = item.split('gene_id')[1].strip().strip('\"')
#             elif 'transcript_id' in item:
#                 transcript_id = item.split('transcript_id')[-1].strip().strip('\"')
#         if not gene_id:
#             logger.error("gene_id not found in '{}'", data.loc[i])
#         if not transcript_id:
#             logger.error("transcript_id not found in '{}'", data.loc[i])
#         # logger.debug("Gene: {}, Transcript: {}", gene_id, transcript_id)
#         data.at[i, "gene_id"] = gene_id
#         data.at[i, "transcript_id"] = transcript_id
#     data = data.drop(labels='attributes', axis=1)
#     logger.debug("Exon data rows: {}", data.shape[0])
#     missing_value_num = data.isnull().sum().sum()
#     if missing_value_num > 0:
#         logger.warning("Total number of missing values: {}", missing_value_num)
#     else:
#         logger.info("No missing values in data")
#     gene_id_missing_value_num = data['gene_id'].isnull().sum()
#     transcript_id_missing_value_num = data['transcript_id'].isnull().sum()
#     if gene_id_missing_value_num > 0:
#         logger.warning("Missing gene_id value number: {}", missing_value_num)
#     if transcript_id_missing_value_num > 0:
#         logger.warning("Missing transcript_id value number: {}", missing_value_num)
#     data['start'] = pd.to_numeric(data['start'], downcast="unsigned")
#     data['end'] = pd.to_numeric(data['end'], downcast="unsigned")
#     return data

def read_exon_data_from_file(infile):
        print ("Reading GTF...")

        omegatic = time.perf_counter()

        
        all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                           'attributes', 'comments']
        
        drop_columns = ['source', 'feature', 'score', 'frame', 'comments']

        data = pd.read_csv(infile, sep='\t', comment='#', header=None, low_memory=False)
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
                
        logger.debug("Exon data rows: {}", newData.shape[0])
        
        missing_value_num = newData.isnull().sum().sum()
        if missing_value_num > 0:
                logger.warning("Total number of missing values: {}", missing_value_num)
        else:
                logger.info("No missing values in data")
        
        gene_id_missing_value_num = newData['gene_id'].isnull().sum()
        
        transcript_id_missing_value_num = newData['transcript_id'].isnull().sum()
        
        if gene_id_missing_value_num > 0:
                logger.warning("Missing gene_id value number: {}", gene_id_missing_value_num)
        if transcript_id_missing_value_num > 0:
                logger.warning("Missing transcript_id value number: {}", transcript_id_missing_value_num)
        
                
        newData['start'] = pd.to_numeric(newData['start'], downcast="unsigned")
        newData['end'] = pd.to_numeric(newData['end'], downcast="unsigned")
        
        toc = time.perf_counter()       
        
        print(f"GTF Read complete,  took {toc-omegatic:0.4f} seconds.")
        return newData


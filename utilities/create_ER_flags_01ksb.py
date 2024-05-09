#!/usr/bin/env python

import argparse
import pandas as pd
import trand.io
import re
import time

def getOptions():
        """
        
        Function to store user input via argparse

        Returns
        -------
        args : ARGPARSE ARGUMENTS
                User input via argparse.

        """
        
        # Parse command line arguments
        parser = argparse.ArgumentParser(description="Identifies groups of transcripts that share all of their exon regions.")
        
        # INPUT
        parser.add_argument(
                "-ea",
                "--event-analysis",
                dest="eaFile",
                required=True,
                help="Location of event analysis ER file from gene mode output"
                )
        
        
        # OUTPUT
        # parser.add_argument("-o",
        #                     "--output-gtf",
        #                     dest="outGTF", 
        #                     required=True,
        #                     help="Path and filename for output GTF")
        
        parser.add_argument("-o",
                            "--output-directory",
                            dest="outDir", 
                            required=True,
                            help="Output path")
        
        args = parser.parse_args()
        return args


def main():
    
    # PROJ="/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/test_gene_based_ERGs"
    PROJ="/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs"
    
    # eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_ea_LOC110191367_noIR.csv"
    # eaFile = "{}/sub_keepIR_event_analysis_er.csv".format(PROJ)
    # eaFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/trand_1GTF_geneMode_fiveSpecies_ujc/fiveSpecies_2_dser1_ujc_event_analysis_er.csv"
    # pw_eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_pw_event_analysis.csv"
    # pw_eaFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/sub_pw_noIR_event_analysis.csv"    
    
    # outGTF = "{}/gene_based_erg_keepIR.gtf".format(PROJ)
    outGTF = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/test_gene_based_ERGs/gene_based_erg_noIR.gtf"
    
    outDir = ""
    
    eaFile = args.eaFile
    outGTF = args.outGTF    
    outDir = args.outDir
    
    eaDf = pd.read_csv(eaFile, low_memory=False)       
    
    tmpDf = eaDf[['er_id','er_transcript_ids','gene_id']].copy()

    # testRowLst = [
    # {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_1', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    # {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_2', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    # {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_3', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    # {'er_id': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268_exon_4', 'er_transcript_ids': '0c10a9919a9cca1fcfd6a40a3e1884251ab221ec924581fdd4ca46a88f951268', 'gene_id': 'LOC110176799'},
    # ]
    # tmpDf = pd.concat([tmpDf, pd.DataFrame(testRowLst)])
    
    tmpDf.columns = ['ER','transcripts','geneID']
    tmpDf['transcripts'] = tmpDf['transcripts'].str.split('|')
    
    tmpDf = tmpDf.explode('transcripts')
    
    perGeneDf = tmpDf.groupby('geneID')
    
    dfLst = []
    wideDfLst = []
    
    alphatic = time.perf_counter()
    for gene, infoDf in perGeneDf:
        
        infoDf['flag_ERPresent'] = 1
        wideDf = pd.pivot_table(infoDf, values='flag_ERPresent', index=['transcripts','geneID'], columns='ER', fill_value=0)
        
        wideDf.columns = sorted(wideDf.columns, key=lambda x: int(x.split("ER")[1]) if 'ER' in x else int(x.split("exon_")[1]))
        
        stackedDf = wideDf.stack().reset_index(name='flag_ERPresent')
        
        dfLst.append(stackedDf)
        wideDfLst.append(wideDf)
    
    outDf = pd.concat(dfLst)
    
    hmm = pd.concat(wideDfLst).fillna(0)
    
    omegatoc = time.perf_counter()
    print(f"Complete! Took {omegatoc-alphatic:0.4f} seconds. Doing next process...")

    outDf.to_csv('{}/test_stacked.csv'.format(outDir),index=False)
    hmm.to_csv('{}/test_wide.csv'.format(outDir),index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
